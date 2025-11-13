// 把你的 IP/LP 實作檔名包含進來（請用你的實際檔名替換這行）
#include "ip_scip.hpp"
#include "sc_params.hpp"
#include <iostream>
#include <chrono>

// ---- 便利函式：統一變數命名（與 IP 的 VarBimap 相容，字串即z是變數 ID）----
static inline std::string vP (const std::string& i, const std::string& j){ return "P[" + i + "," + j + "]"; }
static inline std::string vX (const std::string& i, const std::string& j, const std::string& k){ return "X[" + i + "," + j + "," + k + "]"; }
static inline std::string vY (const std::string& i, const std::string& k, const std::string& l){ return "Y[" + i + "," + k + "," + l + "]"; }
static inline std::string vU (const std::string& i, const std::string& l){ return "U[" + i + "," + l + "]"; }
static inline std::string vW (const std::string& k){ return "W[" + k + "]"; }
static inline std::string vS (const std::string& l){ return "S[" + l + "]"; }

// ---- 核心：依參數建出 IP 模型（目標式 + 限制式）----
IP build_supply_chain_ip(const SCParams& P) {
  const size_t I = P.prod.size();
  const size_t J = P.fac.size();
  const size_t K = P.wh.size();
  const size_t L = P.store.size();

  // ======================
  // 目標式：最大化淨利潤
  // ======================
  std::vector<std::pair<double,std::string>> obj;

  // 銷售收入： + sum_{i,l,k} p_{i,l} * Y_{i,k,l}
  for (size_t i = 0; i < I; ++i)
    for (size_t l = 0; l < L; ++l)
      for (size_t k = 0; k < K; ++k)
        obj.push_back({ + P.price[i][l], vY(P.prod[i], P.wh[k], P.store[l]) });

  // 生產成本： - sum_{i,j} C_{i,j} * P_{i,j}
  for (size_t i = 0; i < I; ++i)
    for (size_t j = 0; j < J; ++j)
      obj.push_back({ - P.prod_cost[i][j], vP(P.prod[i], P.fac[j]) });

  // 物流成本(體積計價)：
  // - sum_{i,j,k} TC1_{j,k} * V_i * X_{i,j,k}
  for (size_t i = 0; i < I; ++i)
    for (size_t j = 0; j < J; ++j)
      for (size_t k = 0; k < K; ++k)
        obj.push_back({ - P.tc1[j][k] * P.V[i], vX(P.prod[i], P.fac[j], P.wh[k]) });

  // - sum_{i,k,l} TC2_{k,l} * V_i * Y_{i,k,l}
  for (size_t i = 0; i < I; ++i)
    for (size_t k = 0; k < K; ++k)
      for (size_t l = 0; l < L; ++l)
        obj.push_back({ - P.tc2[k][l] * P.V[i], vY(P.prod[i], P.wh[k], P.store[l]) });

  // 倉庫/門市固定費： - sum_k R_k W_k  - sum_l SR_l S_l
  for (size_t k = 0; k < K; ++k) obj.push_back({ - P.wh_rent[k],   vW(P.wh[k])   });
  for (size_t l = 0; l < L; ++l) obj.push_back({ - P.store_rent[l], vS(P.store[l]) });

  // 未滿足需求懲罰： - sum_{i,l} M_{i,l} U_{i,l}
  for (size_t i = 0; i < I; ++i)
    for (size_t l = 0; l < L; ++l)
      obj.push_back({ - P.penalty[i][l], vU(P.prod[i], P.store[l]) });

  IP ip("max", obj); // 建立 IP 問題（純整數：所有變數皆為非負整數）

  // =================================
  // 限制式群組
  // =================================

  // (1) 工廠工時產能： sum_i T_{i,j} P_{i,j} <= Cap_j
  for (size_t j = 0; j < J; ++j) {
    std::vector<std::pair<double,std::string>> terms;
    for (size_t i = 0; i < I; ++i)
      terms.push_back({ P.prod_time[i][j], vP(P.prod[i], P.fac[j]) });
    ip.addConstraint(terms, "<=", P.cap[j]);
  }

  // (2) 生產 = 出廠： P_{i,j} - sum_k X_{i,j,k} = 0
  for (size_t i = 0; i < I; ++i)
    for (size_t j = 0; j < J; ++j) {
      std::vector<std::pair<double,std::string>> terms;
      terms.push_back({ +1.0, vP(P.prod[i], P.fac[j]) });
      for (size_t k = 0; k < K; ++k)
        terms.push_back({ -1.0, vX(P.prod[i], P.fac[j], P.wh[k]) });
      ip.addConstraint(terms, "=", 0.0);
    }

  // (3) 倉庫流量守恆： sum_j X_{i,j,k} - sum_l Y_{i,k,l} = 0
  for (size_t i = 0; i < I; ++i)
    for (size_t k = 0; k < K; ++k) {
      std::vector<std::pair<double,std::string>> terms;
      for (size_t j = 0; j < J; ++j)
        terms.push_back({ +1.0, vX(P.prod[i], P.fac[j], P.wh[k]) });
      for (size_t l = 0; l < L; ++l)
        terms.push_back({ -1.0, vY(P.prod[i], P.wh[k], P.store[l]) });
      ip.addConstraint(terms, "=", 0.0);
    }

  // (4) 倉庫吞吐容量（體積）與啟用邏輯：
  // sum_i V_i * (sum_j X_{i,j,k}) - SCap_k * W_k <= 0
  for (size_t k = 0; k < K; ++k) {
    std::vector<std::pair<double,std::string>> terms;
    for (size_t i = 0; i < I; ++i)
      for (size_t j = 0; j < J; ++j)
        terms.push_back({ + P.V[i], vX(P.prod[i], P.fac[j], P.wh[k]) });
    terms.push_back({ - P.wh_cap[k], vW(P.wh[k]) });
    ip.addConstraint(terms, "<=", 0.0);
  }

  // (5) 市場需求平衡（含未滿足）：
  // sum_k Y_{i,k,l} + U_{i,l} = D_{i,l}
  for (size_t i = 0; i < I; ++i)
    for (size_t l = 0; l < L; ++l) {
      std::vector<std::pair<double,std::string>> terms;
      for (size_t k = 0; k < K; ++k)
        terms.push_back({ +1.0, vY(P.prod[i], P.wh[k], P.store[l]) });
      terms.push_back({ +1.0, vU(P.prod[i], P.store[l]) });
      ip.addConstraint(terms, "=", P.demand[i][l]);
    }

  // (6) 未滿足需求上界（強化）： U_{i,l} <= D_{i,l}
  for (size_t i = 0; i < I; ++i)
    for (size_t l = 0; l < L; ++l) {
      std::vector<std::pair<double,std::string>> terms = { {+1.0, vU(P.prod[i], P.store[l])} };
      ip.addConstraint(terms, "<=", P.demand[i][l]);
    }

  // (7) 門市啟用邏輯（Big-M = D_{i,l}）： sum_k Y_{i,k,l} - D_{i,l} * S_l <= 0
  for (size_t i = 0; i < I; ++i)
    for (size_t l = 0; l < L; ++l) {
      std::vector<std::pair<double,std::string>> terms;
      for (size_t k = 0; k < K; ++k)
        terms.push_back({ +1.0, vY(P.prod[i], P.wh[k], P.store[l]) });
      terms.push_back({ - P.demand[i][l], vS(P.store[l]) });
      ip.addConstraint(terms, "<=", 0.0);
    }

  // (8) 二元邏輯：W_k, S_l ∈ {0,1}  ——> 0 ≤ var ≤ 1（IP 預設非負整數，故加上 ≤1 即為二元）
  for (size_t k = 0; k < K; ++k) {
    std::vector<std::pair<double,std::string>> terms = { {+1.0, vW(P.wh[k])} };
    ip.addConstraint(terms, "<=", 1.0);
  }
  for (size_t l = 0; l < L; ++l) {
    std::vector<std::pair<double,std::string>> terms = { {+1.0, vS(P.store[l])} };
    ip.addConstraint(terms, "<=", 1.0);
  }

  return ip;
}

int main() {
    // 設定參數
    SCGenCfg cfg;
    cfg.I = 5;   // 5 種產品
    cfg.J = 3;   // 3 個工廠
    cfg.K = 2;   // 2 個倉庫
    cfg.L = 4;   // 4 個門市
    
    SCParams params = make_sc_params(cfg);
    
    std::cout << "Building model..." << std::endl;
    std::cout << "Products: " << cfg.I << ", Factories: " << cfg.J 
              << ", Warehouses: " << cfg.K << ", Stores: " << cfg.L << std::endl;
    
    // 建立模型
    IP model = build_supply_chain_ip(params);
    
    // 設定 SCIP 參數
    model.setTimeLimit(300.0);  // 5 分鐘時間限制
    model.setVerbosity(4);      // 詳細輸出
    
    // 求解
    std::cout << "\nSolving with SCIP..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    bool success = model.solve();
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // 輸出結果
    std::cout << "\n==================== Results ====================" << std::endl;
    std::cout << "Status: " << model.getStatus() << std::endl;
    std::cout << "Solve time: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    if (success) {
        std::cout << "Objective value: " << model.getObjValue() << std::endl;
        std::cout << "\nSolution details:" << std::endl;
        model.printSolution();
    } else {
        std::cout << "No feasible solution found!" << std::endl;
    }
    
    return 0;
}