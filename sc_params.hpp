#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>      // floor, round

// ===================
// 參數資料結構（沿用你的定義，以 double 存整數值）
// ===================
struct SCParams {
  std::vector<std::string> prod;   // I: 產品
  std::vector<std::string> fac;    // J: 工廠
  std::vector<std::string> wh;     // K: 倉庫
  std::vector<std::string> store;  // L: 門市

  std::vector<double> V; // |I| 產品體積（立方公尺/件）

  std::vector<std::vector<double>> price;   // |I| x |L| 售價
  std::vector<std::vector<double>> demand;  // |I| x |L| 需求上限
  std::vector<std::vector<double>> penalty; // |I| x |L| 未滿足懲罰

  std::vector<std::vector<double>> prod_cost; // |I| x |J| 生產成本
  std::vector<std::vector<double>> prod_time; // |I| x |J| 單位工時
  std::vector<double> cap; // |J| 工廠工時上限

  std::vector<double> wh_rent; // |K| 倉庫固定費
  std::vector<double> wh_cap;  // |K| 倉庫吞吐容量（以體積計）

  std::vector<double> store_rent; // |L| 門市固定費

  std::vector<std::vector<double>> tc1; // |J| x |K| 工廠→倉庫（單位體積）
  std::vector<std::vector<double>> tc2; // |K| x |L| 倉庫→門市（單位體積）
};

// -------------------
// 命名工具
// -------------------
inline std::vector<std::string> make_product_names(int n) {
  std::vector<std::string> v;
  v.reserve(n);
  for (int i = 0; i < n; ++i) {
    char base = 'A' + (i % 26);
    int round = i / 26;
    if (round == 0) v.emplace_back(std::string(1, base));
    else v.emplace_back(std::string(1, base) + std::to_string(round + 1));
  }
  return v;
}
inline std::vector<std::string> make_seq_names(const char* prefix, int n) {
  std::vector<std::string> v;
  v.reserve(n);
  for (int i = 0; i < n; ++i) v.emplace_back(std::string(prefix) + std::to_string(i+1));
  return v;
}

// ===================
// 參數生成「可調」超參數
//  - 設計原則：數值偏小但合理、整數化、且保證單件正利潤
// ===================
struct SCGenCfg {
  // 尺寸
  int I = 3, J = 2, K = 1, L = 2;

  // 體積 V_i = vol_start + i * vol_step
  int vol_start = 1;
  int vol_step  = 1;

  // 單位工時 T_{i,j} = time_base + i + (j%2)*time_parity_bonus
  int time_base = 1;
  int time_parity_bonus = 1;

  // 生產成本 base_cost_i = cost_base + cost_step * i
  // 各工廠再乘上 (100 + shift)/100，其中 shift 線性分佈於 [-cost_grad_pct, +cost_grad_pct]
  int cost_base = 200;
  int cost_step = 100;
  int cost_grad_pct = 8; // 0~100, 代表 ±%

  // 需求 D_{i,l} = demand_base + demand_i_step*i + demand_l_step*(l%4)
  int demand_base    = 20;
  int demand_i_step  = 5;
  int demand_l_step  = 3;

  // 運費（以體積計）：TC1_{j,k} = tc1_base + tc_step*((j%3)+(k%4))
  //                     TC2_{k,l} = tc2_base + tc_step*((k%4)+(l%4))
  int tc1_base = 8;
  int tc2_base = 9;
  int tc_step  = 2;

  // 售價 price[i][l] = minProd_i + V_i * minShipPerVol_l + margin_i
  // margin_i = max( floor(margin_frac * minProd_i), margin_floor_base + margin_floor_step*i )
  double margin_frac = 0.25; // 25% of min production cost
  int    margin_floor_base = 20;
  int    margin_floor_step = 5;

  // 未滿足懲罰：penalty = floor(penalty_frac * price)
  double penalty_frac = 0.6;

  // 工廠工時上限：Cap_j ≈ cap_util × (總需求工時 / J) + cap_buffer
  double cap_util   = 0.7; // 70% 目標稼動
  int    cap_buffer = 50;  // 小緩衝

  // 倉庫吞吐容量：wh_cap ≈ wh_capacity_share × (總需求體積 / K)，至少 1
  double wh_capacity_share = 0.5;

  // 固定費（小）：避免壓過利潤
  int wh_rent_base   = 2000;
  int wh_rent_step   = 200;
  int store_rent_base= 6000;
  int store_rent_step= 500;
};

// ===================
// 主要生成器：依 SCGenCfg 產生 SCParams（整數化 + 正利潤）
// ===================
inline SCParams make_sc_params(const SCGenCfg& C) {
  SCParams P;

  // ---- 集合命名 ----
  P.prod  = make_product_names(C.I);
  P.fac   = make_seq_names("F", C.J);
  P.wh    = make_seq_names("W", C.K);
  P.store = make_seq_names("S", C.L);

  // ---- 體積 V_i（整數）----
  P.V.assign(C.I, 0.0);
  for (int i = 0; i < C.I; ++i) {
    int v = std::max(1, C.vol_start + C.vol_step * i);
    P.V[i] = static_cast<double>(v);
  }

  // ---- 單位工時 T_{i,j}（整數）----
  P.prod_time.assign(C.I, std::vector<double>(C.J, 0.0));
  for (int i = 0; i < C.I; ++i)
    for (int j = 0; j < C.J; ++j) {
      int t = C.time_base + i + (j % 2) * C.time_parity_bonus;
      P.prod_time[i][j] = static_cast<double>(std::max(1, t));
    }

  // ---- 生產成本 C_{i,j}（整數）----
  // base_cost_i = cost_base + cost_step * i；工廠依位置給 ±cost_grad_pct 線性梯度
  P.prod_cost.assign(C.I, std::vector<double>(C.J, 0.0));
  for (int i = 0; i < C.I; ++i) {
    int base_cost_i = std::max(1, C.cost_base + C.cost_step * i);
    for (int j = 0; j < C.J; ++j) {
      int shift = 0;
      if (C.J > 1) {
        // 線性由 -grad -> +grad
        shift = (j * (2 * C.cost_grad_pct)) / (C.J - 1) - C.cost_grad_pct;
      }
      int cost_ij = base_cost_i * (100 + shift) / 100;
      P.prod_cost[i][j] = static_cast<double>(std::max(1, cost_ij));
    }
  }

  // ---- 需求 D_{i,l}（整數）----
  P.demand.assign(C.I, std::vector<double>(C.L, 0.0));
  for (int i = 0; i < C.I; ++i)
    for (int l = 0; l < C.L; ++l) {
      int d = C.demand_base + C.demand_i_step * i + C.demand_l_step * (l % 4);
      P.demand[i][l] = static_cast<double>(std::max(0, d));
    }

  // ---- 運費（整數、以體積計）----
  P.tc1.assign(C.J, std::vector<double>(C.K, 0.0));
  for (int j = 0; j < C.J; ++j)
    for (int k = 0; k < C.K; ++k) {
      int v = C.tc1_base + C.tc_step * ((j % 3) + (k % 4));
      P.tc1[j][k] = static_cast<double>(std::max(0, v));
    }

  P.tc2.assign(C.K, std::vector<double>(C.L, 0.0));
  for (int k = 0; k < C.K; ++k)
    for (int l = 0; l < C.L; ++l) {
      int v = C.tc2_base + C.tc_step * ((k % 4) + (l % 4));
      P.tc2[k][l] = static_cast<double>(std::max(0, v));
    }

  // ---- 計算 min 生產成本 & 最便宜單位體積運價（保證正利潤）----
  std::vector<int> minProd(C.I, 0);
  for (int i = 0; i < C.I; ++i) {
    int mn = std::numeric_limits<int>::max();
    for (int j = 0; j < C.J; ++j) mn = std::min(mn, static_cast<int>(P.prod_cost[i][j]));
    minProd[i] = mn;
  }
  std::vector<int> minShipPerVol(C.L, 0);
  for (int l = 0; l < C.L; ++l) {
    int best = std::numeric_limits<int>::max();
    for (int k = 0; k < C.K; ++k) {
      int bestF = std::numeric_limits<int>::max();
      for (int j = 0; j < C.J; ++j)
        bestF = std::min(bestF, static_cast<int>(P.tc1[j][k]));
      int total = bestF + static_cast<int>(P.tc2[k][l]);
      best = std::min(best, total);
    }
    minShipPerVol[l] = best; // 每立方公尺的最便宜路徑成本
  }

  // ---- 售價（整數，確保單件正利潤）----
  P.price.assign(C.I, std::vector<double>(C.L, 0.0));
  for (int i = 0; i < C.I; ++i) {
    int base_margin = static_cast<int>(std::floor(minProd[i] * C.margin_frac));
    base_margin = std::max(base_margin, C.margin_floor_base + C.margin_floor_step * i);
    base_margin = std::max(1, base_margin);
    for (int l = 0; l < C.L; ++l) {
      int ship = static_cast<int>(P.V[i]) * std::max(0, minShipPerVol[l]);
      int price_il = minProd[i] + ship + base_margin;
      // 下限：至少比 (minProd + ship) 多 1，避免 margin=0
      price_il = std::max(price_il, minProd[i] + ship + 1);
      P.price[i][l] = static_cast<double>(price_il);
    }
  }

  // ---- 未滿足懲罰（整數）----
  P.penalty.assign(C.I, std::vector<double>(C.L, 0.0));
  for (int i = 0; i < C.I; ++i)
    for (int l = 0; l < C.L; ++l) {
      int pen = static_cast<int>(std::floor(P.price[i][l] * C.penalty_frac));
      P.penalty[i][l] = static_cast<double>(std::max(0, pen));
    }

  // ---- 工廠工時上限 Cap_j（整數）----
  std::vector<int> sumD(C.I, 0);
  for (int i = 0; i < C.I; ++i) {
    int s = 0; for (int l = 0; l < C.L; ++l) s += static_cast<int>(P.demand[i][l]);
    sumD[i] = s;
  }
  P.cap.assign(C.J, 0.0);
  for (int j = 0; j < C.J; ++j) {
    long long hours = 0;
    for (int i = 0; i < C.I; ++i)
      hours += 1LL * sumD[i] * static_cast<int>(P.prod_time[i][j]);
    long long capj = static_cast<long long>(std::floor((hours / std::max(1, C.J)) * C.cap_util))
                   + C.cap_buffer;
    P.cap[j] = static_cast<double>(std::max<long long>(1, capj));
  }

  // ---- 倉庫吞吐容量 wh_cap（整數體積）----
  long long totalVol = 0;
  for (int i = 0; i < C.I; ++i)
    totalVol += 1LL * sumD[i] * static_cast<int>(P.V[i]);
  P.wh_cap.assign(C.K, 0.0);
  for (int k = 0; k < C.K; ++k) {
    long long capk = static_cast<long long>(std::floor((totalVol * C.wh_capacity_share) / std::max(1, C.K)));
    P.wh_cap[k] = static_cast<double>(std::max<long long>(1, capk));
  }

  // ---- 固定費（整數，小）----
  P.wh_rent.assign(C.K, 0.0);
  for (int k = 0; k < C.K; ++k)
    P.wh_rent[k] = static_cast<double>(C.wh_rent_base + C.wh_rent_step * (k + 1));

  P.store_rent.assign(C.L, 0.0);
  for (int l = 0; l < C.L; ++l)
    P.store_rent[l] = static_cast<double>(C.store_rent_base + C.store_rent_step * (l + 1));

  return P;
}

// -------------------
// 便利介面：保留舊版 API（預設尺寸）
// -------------------
inline SCParams default_sc_params(int I = 2, int J = 2, int K = 1, int L = 2) {
  SCGenCfg cfg;
  cfg.I = I; cfg.J = J; cfg.K = K; cfg.L = L;
  return make_sc_params(cfg);
}
