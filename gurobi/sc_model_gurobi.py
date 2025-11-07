# -*- coding: utf-8 -*-
"""
三階層供應鏈（工廠→倉庫→門市），單期、純整數變數的 MIP
- 目標：最大化淨利潤
- 限制：工時產能、流量守恆、倉庫容量(體積)且需啟用、需求平衡(含未滿足)、門市啟用邏輯
可用參數：
  --relax        ：把整數流量放寬為連續（W/S 也放寬為 [0,1]），做 LP 對比
  --timelimit T  ：秒
  --focus k      ：MIPFocus（0 預設、1 可行解、2 下界、3 最佳性）
"""
import argparse
import gurobipy as gp
from gurobipy import GRB
from sc_params import get_default_params

def build_and_solve(relax: bool=False, timelimit: float=None, mipfocus: int=None):
    P = get_default_params()
    I, J, K, L = P["I"], P["J"], P["K"], P["L"]
    V = P["V"]; price=P["price"]; demand=P["demand"]; penalty=P["penalty"]
    prod_cost=P["prod_cost"]; prod_time=P["prod_time"]; cap=P["cap"]
    wh_rent=P["wh_rent"]; wh_cap=P["wh_cap"]; store_rent=P["store_rent"]
    tc1=P["tc1"]; tc2=P["tc2"]

    m = gp.Model("SupplyChain_IP")

    # ============= 變數 =============
    # 依 --relax 控制型別
    flow_vtype = GRB.CONTINUOUS if relax else GRB.INTEGER
    bin_vtype  = GRB.CONTINUOUS if relax else GRB.BINARY
    bin_ub     = 1.0 if relax else 1.0  # 放寬時改為 0..1

    # P[i,j] 產量
    Pvar = m.addVars(I, J, name="P", vtype=flow_vtype, lb=0)

    # X[i,j,k] 工廠→倉庫 出貨
    Xvar = m.addVars(I, J, K, name="X", vtype=flow_vtype, lb=0)

    # Y[i,k,l] 倉庫→門市 出貨
    Yvar = m.addVars(I, K, L, name="Y", vtype=flow_vtype, lb=0)

    # U[i,l] 未滿足需求（上界 = D[i,l]）
    Uvar = m.addVars(I, L, name="U", vtype=flow_vtype, lb=0,
                     ub={(i,l): demand[i][l] for i in I for l in L})

    # W[k] 倉庫啟用、S[l] 門市開店
    Wvar = m.addVars(K, name="W", vtype=bin_vtype, lb=0, ub=bin_ub)
    Svar = m.addVars(L, name="S", vtype=bin_vtype, lb=0, ub=bin_ub)

    # ============= 目標式 =============
    obj = gp.LinExpr()

    # 銷售收入： + sum_{i,l,k} p_{i,l} * Y_{i,k,l}
    obj += gp.quicksum(price[i][l] * Yvar[i,k,l] for i in I for l in L for k in K)

    # 生產成本： - sum_{i,j} C_{i,j} * P_{i,j}
    obj += gp.quicksum((-prod_cost[i][j]) * Pvar[i,j] for i in I for j in J)

    # 物流成本（體積計價）：
    # - sum_{i,j,k} TC1_{j,k} * V_i * X_{i,j,k}
    obj += gp.quicksum((-tc1[j][k] * V[i]) * Xvar[i,j,k] for i in I for j in J for k in K)

    # - sum_{i,k,l} TC2_{k,l} * V_i * Y_{i,k,l}
    obj += gp.quicksum((-tc2[k][l] * V[i]) * Yvar[i,k,l] for i in I for k in K for l in L)

    # 固定費： - sum_k R_k W_k - sum_l SR_l S_l
    obj += gp.quicksum((-wh_rent[k]) * Wvar[k] for k in K)
    obj += gp.quicksum((-store_rent[l]) * Svar[l] for l in L)

    # 未滿足懲罰： - sum_{i,l} M_{i,l} U_{i,l}
    obj += gp.quicksum((-penalty[i][l]) * Uvar[i,l] for i in I for l in L)

    m.setObjective(obj, GRB.MAXIMIZE)

    # ============= 限制式 =============

    # (1) 工廠工時產能： sum_i T_{i,j} P_{i,j} <= Cap_j
    m.addConstrs(
        (gp.quicksum(prod_time[i][j]*Pvar[i,j] for i in I) <= cap[j]
         for j in J),
        name="FactoryTime"
    )

    # (2) 生產 = 出廠： P_{i,j} = sum_k X_{i,j,k}
    m.addConstrs(
        (Pvar[i,j] == gp.quicksum(Xvar[i,j,k] for k in K)
         for i in I for j in J),
        name="ProdBalance"
    )

    # (3) 倉庫流量守恆： sum_j X_{i,j,k} = sum_l Y_{i,k,l}
    m.addConstrs(
        (gp.quicksum(Xvar[i,j,k] for j in J) == gp.quicksum(Yvar[i,k,l] for l in L)
         for i in I for k in K),
        name="WhBalance"
    )

    # (4) 倉庫吞吐容量(體積)與啟用： sum_i V_i * sum_j X_{i,j,k} <= SCap_k * W_k
    m.addConstrs(
        (gp.quicksum(V[i]*gp.quicksum(Xvar[i,j,k] for j in J) for i in I)
         <= wh_cap[k] * Wvar[k]
         for k in K),
        name="WhCapActivate"
    )

    # (5) 需求平衡： sum_k Y_{i,k,l} + U_{i,l} = D_{i,l}
    m.addConstrs(
        (gp.quicksum(Yvar[i,k,l] for k in K) + Uvar[i,l] == demand[i][l]
         for i in I for l in L),
        name="DemandBalance"
    )

    # (6) 門市啟用邏輯（Big-M = D_{i,l}）： sum_k Y_{i,k,l} <= D_{i,l} * S_l
    m.addConstrs(
        (gp.quicksum(Yvar[i,k,l] for k in K) <= demand[i][l]*Svar[l]
         for i in I for l in L),
        name="StoreActivate"
    )

    # 一些解器參數（可依需求調整）
    if timelimit is not None:
        m.Params.TimeLimit = float(timelimit)
    if mipfocus is not None:
        m.Params.MIPFocus = int(mipfocus)
    # m.Params.OutputFlag = 1  # 需要時顯示更詳細日誌

    # 求解
    m.optimize()

    # ============= 輸出重點結果 =============
    print("\n=== Solve Summary ===")
    status_map = {
        GRB.OPTIMAL: "OPTIMAL",
        GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INTERRUPTED: "INTERRUPTED",
        GRB.INFEASIBLE: "INFEASIBLE",
        GRB.UNBOUNDED: "UNBOUNDED",
        GRB.INF_OR_UNBD: "INF_OR_UNBD",
    }
    print("Status :", status_map.get(m.Status, f"CODE_{m.Status}"))
    if m.SolCount > 0:
        print("ObjVal :", m.ObjVal)
        if m.Status != GRB.INFEASIBLE and m.Status != GRB.UNBOUNDED and m.Status != GRB.INF_OR_UNBD:
            try:
                print("MIPGap :", m.MIPGap)
            except Exception:
                pass

        # 啟用的倉庫/門市
        open_wh = [k for k in K if Wvar[k].X > 0.5]
        open_st = [l for l in L if Svar[l].X > 0.5]
        print("Open Warehouses:", open_wh)
        print("Open Stores    :", open_st)

        # 顯示非零決策（節省輸出，> 1e-6 才列）
        thr = 1e-6
        def nz_items(vardict, fmt):
            for idx, var in vardict.items():
                if var.X > thr:
                    print(fmt.format(idx=idx, val=var.X))

        print("\n-- Production P[i,j] --")
        for i in I:
            for j in J:
                val = Pvar[i,j].X
                if val > thr:
                    print(f"P[{i},{j}] = {val:.4f}")

        print("\n-- Flow X[i,j,k] --")
        for i in I:
            for j in J:
                for k in K:
                    val = Xvar[i,j,k].X
                    if val > thr:
                        print(f"X[{i},{j},{k}] = {val:.4f}")

        print("\n-- Flow Y[i,k,l] --")
        for i in I:
            for k in K:
                for l in L:
                    val = Yvar[i,k,l].X
                    if val > thr:
                        print(f"Y[{i},{k},{l}] = {val:.4f}")

        print("\n-- Unmet U[i,l] --")
        for i in I:
            for l in L:
                val = Uvar[i,l].X
                if val > thr:
                    print(f"U[{i},{l}] = {val:.4f}")

    return m

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--relax", action="store_true",
                    help="將整數模型放寬為 LP 鬆弛（W/S 變為 0..1 連續），做對比測試")
    ap.add_argument("--timelimit", type=float, default=None, help="時間上限（秒）")
    ap.add_argument("--focus", type=int, default=None, help="MIPFocus ∈ {0,1,2,3}")
    args = ap.parse_args()
    build_and_solve(relax=args.relax, timelimit=args.timelimit, mipfocus=args.focus)
