# -*- coding: utf-8 -*-
"""
自動生成參數（Python 版）
- 介面乾淨：透過函式參數調整超參數（與 C++ SCGenCfg 一致的語意）
- 數值偏小、整數化、且保證每件正利潤（price > minProd + V*minShipPerVol）
"""

from typing import Dict, List


def _make_product_names(n: int) -> List[str]:
    v = []
    for i in range(n):
        base = chr(ord('A') + (i % 26))
        round_ = i // 26
        v.append(base if round_ == 0 else f"{base}{round_+1}")
    return v


def _make_seq(prefix: str, n: int) -> List[str]:
    return [f"{prefix}{i+1}" for i in range(n)]


def gen_params(
    I: int = 3, J: int = 2, K: int = 1, L: int = 2,
    # 體積
    vol_start: int = 1, vol_step: int = 1,
    # 工時
    time_base: int = 1, time_parity_bonus: int = 1,
    # 生產成本
    cost_base: int = 200, cost_step: int = 100, cost_grad_pct: int = 8,
    # 需求
    demand_base: int = 20, demand_i_step: int = 5, demand_l_step: int = 3,
    # 運費（以體積計）
    tc1_base: int = 8, tc2_base: int = 9, tc_step: int = 2,
    # 價格與懲罰
    margin_frac: float = 0.25, margin_floor_base: int = 20, margin_floor_step: int = 5,
    penalty_frac: float = 0.6,
    # 產能與倉容
    cap_util: float = 0.7, cap_buffer: int = 50, wh_capacity_share: float = 0.5,
    # 固定費（小）
    wh_rent_base: int = 2000, wh_rent_step: int = 200,
    store_rent_base: int = 6000, store_rent_step: int = 500,
):
    # === 集合 ===
    I_names: List[str] = _make_product_names(I)
    J_names: List[str] = _make_seq("F", J)
    K_names: List[str] = _make_seq("W", K)
    L_names: List[str] = _make_seq("S", L)

    # === 體積（整數）V_i = vol_start + vol_step*i，至少 1 ===
    V: Dict[str, int] = {}
    for i, pi in enumerate(I_names):
        V[pi] = max(1, vol_start + vol_step * i)

    # === 單位工時（整數）T_{i,j} = time_base + i + (j%2)*time_parity_bonus，至少 1 ===
    prod_time: Dict[str, Dict[str, int]] = {}
    for i, pi in enumerate(I_names):
        row = {}
        for j, fj in enumerate(J_names):
            t = time_base + i + (j % 2) * time_parity_bonus
            row[fj] = max(1, t)
        prod_time[pi] = row

    # === 生產成本（整數）base_cost_i = cost_base + cost_step*i，工廠 ±cost_grad_pct 線性梯度 ===
    prod_cost: Dict[str, Dict[str, int]] = {}
    for i, pi in enumerate(I_names):
        base = max(1, cost_base + cost_step * i)
        row = {}
        for j, fj in enumerate(J_names):
            if J > 1:
                shift = (j * (2 * cost_grad_pct)) // (J - 1) - cost_grad_pct
            else:
                shift = 0
            row[fj] = max(1, base * (100 + shift) // 100)
        prod_cost[pi] = row

    # === 需求（整數）D_{i,l} = demand_base + demand_i_step*i + demand_l_step*(l%4) ===
    demand: Dict[str, Dict[str, int]] = {}
    for i, pi in enumerate(I_names):
        row = {}
        for l, sl in enumerate(L_names):
            d = demand_base + demand_i_step * i + demand_l_step * (l % 4)
            row[sl] = max(0, d)
        demand[pi] = row

    # === 運費（整數、以體積計）===
    tc1: Dict[str, Dict[str, int]] = {}
    for j, fj in enumerate(J_names):
        row = {}
        for k, wk in enumerate(K_names):
            v = tc1_base + tc_step * ((j % 3) + (k % 4))
            row[wk] = max(0, v)
        tc1[fj] = row

    tc2: Dict[str, Dict[str, int]] = {}
    for k, wk in enumerate(K_names):
        row = {}
        for l, sl in enumerate(L_names):
            v = tc2_base + tc_step * ((k % 4) + (l % 4))
            row[sl] = max(0, v)
        tc2[wk] = row

    # === min 生產成本 & 每門市最便宜單位體積運價 ===
    min_prod: Dict[str, int] = {pi: min(prod_cost[pi][fj] for fj in J_names) for pi in I_names}

    min_ship_per_vol: Dict[str, int] = {}
    for sl in L_names:
        best = 10**9
        for wk in K_names:
            bestF = min(tc1[fj][wk] for fj in J_names)
            best = min(best, bestF + tc2[wk][sl])
        min_ship_per_vol[sl] = best if best < 10**9 else 0  # K=0 不會發生，保險

    # === 售價（整數、保證正利潤） price = minProd + V * minShipPerVol + margin_i ===
    price: Dict[str, Dict[str, int]] = {}
    for i, pi in enumerate(I_names):
        base_margin = max(int(min_prod[pi] * margin_frac), margin_floor_base + margin_floor_step * i, 1)
        row = {}
        for sl in L_names:
            ship = V[pi] * max(0, min_ship_per_vol[sl])
            p = min_prod[pi] + ship + base_margin
            p = max(p, min_prod[pi] + ship + 1)  # 強制至少 +1，確保 margin>0
            row[sl] = p
        price[pi] = row

    # === 未滿足懲罰（整數） penalty = floor(penalty_frac * price) ===
    penalty: Dict[str, Dict[str, int]] = {
        pi: {sl: int(price[pi][sl] * penalty_frac) for sl in L_names} for pi in I_names
    }

    # === 工廠工時上限（整數） Cap ≈ cap_util × (總需求工時 / J) + cap_buffer ===
    sumD: Dict[str, int] = {pi: sum(demand[pi][sl] for sl in L_names) for pi in I_names}
    cap: Dict[str, int] = {}
    for fj in J_names:
        hours = 0
        for pi in I_names:
            hours += sumD[pi] * prod_time[pi][fj]
        cap[fj] = max(1, int((hours / max(1, J)) * cap_util) + cap_buffer)

    # === 倉庫吞吐容量（整數體積）≈ wh_capacity_share × (總需求體積 / K) ===
    total_vol = 0
    for pi in I_names:
        total_vol += sumD[pi] * V[pi]
    wh_cap: Dict[str, int] = {}
    for wk in K_names:
        capk = int((total_vol * wh_capacity_share) / max(1, K))
        wh_cap[wk] = max(1, capk)

    # === 固定費（小、整數） ===
    wh_rent: Dict[str, int]   = {wk: wh_rent_base   + wh_rent_step   * (i+1) for i, wk in enumerate(K_names)}
    store_rent: Dict[str, int]= {sl: store_rent_base+ store_rent_step* (i+1) for i, sl in enumerate(L_names)}

    return dict(
        I=I_names, J=J_names, K=K_names, L=L_names,
        V=V, price=price, demand=demand, penalty=penalty,
        prod_cost=prod_cost, prod_time=prod_time, cap=cap,
        wh_rent=wh_rent, wh_cap=wh_cap, store_rent=store_rent,
        tc1=tc1, tc2=tc2
    )


def get_default_params(
    I: int = 2,
    J: int = 2,
    K: int = 1,
    L: int = 2,
    **overrides
):
    """
    取得「預設參數集」——語意等同於 C++ 的 default_sc_params(I,J,K,L)。

    參數
    ----
    I, J, K, L : int
        產品/工廠/倉庫/門市數量。
    **overrides :
        若需要，還可同時覆蓋 gen_params() 內其他超參數，例如：
        margin_frac=0.3, tc1_base=6, tc2_base=7, demand_base=15, ...

    回傳
    ----
    dict
        與 gen_params() 相同結構：
        {
          "I": [...], "J": [...], "K": [...], "L": [...],
          "V": {...}, "price": {...}, "demand": {...}, "penalty": {...},
          "prod_cost": {...}, "prod_time": {...}, "cap": {...},
          "wh_rent": {...}, "wh_cap": {...}, "store_rent": {...},
          "tc1": {...}, "tc2": {...}
        }
    """
    return gen_params(I=I, J=J, K=K, L=L, **overrides)


# 方便的別名（若你想用和 C++ 一樣的函式名）
default_sc_params = get_default_params
