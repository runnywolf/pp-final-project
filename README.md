# 2025 NYCU Parallel Programming Final Project
Author:
[@homer40214](https://github.com/homer40214)
[@Jerry0223](https://github.com/Jerry0223)
[@runnywolf](https://github.com/runnywolf)

## How to run
Need `c++17`

Login `<username>@hpclogin[01~03].cs.nycu.edu.tw`

Then
```sh
git clone https://github.com/runnywolf/pp-final-project.git
cd pp-final-project
```

Run
```sh
make && ./main.out
```





### 新增/調整的檔案結構

```
pp-final-project/
├─ sc_params.hpp            # C++ 參數（多產品/工廠/倉庫/門市）
├─ sc_model.cpp             # C++ 限制式/目標式建模（餵給 IP 解器）
├─ gurobi/
│  ├─ sc_params.py          # Python 參數（同上，但給 gurobipy 用）
│  └─ sc_model_gurobi.py    # Python 建模與求解（限制式 + 目標式 + 執行）
├─ main.cpp / ip_solver.hpp # 你的 IP/LP 解器與主程式
└─ ...
```

* **C++ 路徑**

  * `sc_params.hpp`：集中所有參數（集合、售價、需求、懲罰、生產成本/工時、容量、固定費、運費…）。
  * `sc_model.cpp`：把上面參數轉成**目標式 + 限制式**，建出 `IP` 物件（純整數）給你的解器。

* **Gurobi（Python）路徑**

  * `gurobi/sc_params.py`：與 `sc_params.hpp` 對齊的參數檔。
  * `gurobi/sc_model_gurobi.py`：用 `gurobipy` 建模與求解（支援 `--relax` 做 LP 鬆弛對比）。

---
## 用 Gurobi 驗證（Python 版）

本專案同時提供 **C++ 自製 IP 解器** 與 **Gurobi（gurobipy）** 兩套模型，方便對解驗證。

### 0) 學生帳號註冊（已完成可跳過）

1. 到 [https://www.gurobi.com](https://www.gurobi.com) 註冊帳號，申請 **Free Academic**。
2. 建議使用 **Academic WLS（Web License）**，免安裝 `grbgetkey`，較不受 HostID 影響。

---

### 1)（建議）WLS 授權設定（免 `grbgetkey`）

> **WLS**：雲端授權，適合學校機器/WSL/容器環境。求解時需可連網。

1. 到 **Gurobi License Portal**（Web License Manager）建立 **WLS API Key**，下載 `gurobi.lic`（內容包含 `WLSACCESSID` / `WLSSECRET` / `LICENSEID` 三行）。
2. 把 `gurobi.lic` 放到登入帳號的家目錄：

   ```sh
   mkdir -p ~/.gurobi
   mv /path/to/gurobi.lic ~/.gurobi/gurobi.lic
   ```

   > 若你之前有設定 `GRB_LICENSE_FILE` 指到別處，建議先 `unset GRB_LICENSE_FILE`，讓 Gurobi 走預設路徑。

---

### 2) 在 hpclogin 建 Python 環境並安裝 `gurobipy`

> **提示**：不同 hpclogin 節點 Python 版本可能略有差異，建議自建虛擬環境。

```sh
# 建議在專案根目錄
python3 -m venv .venv-gurobi
source .venv-gurobi/bin/activate            # (Windows/WSL 改用: .\.venv-gurobi\Scripts\activate)
python -m pip install --upgrade pip
python -m pip install gurobipy
```

（選擇性）驗證安裝與授權：

```sh
python - << 'PY'
import gurobipy as gp
m = gp.Model("check")
m.optimize()
print("gurobipy =", gp.__version__, "| solved with status", m.Status)
PY
```

---

### 3) 執行 Gurobi 模型（Python）

在專案根目錄（已放好 `gurobi.lic` 且啟用虛擬環境）：

```sh
# 整數模型（純整數）
python gurobi/sc_model_gurobi.py

# LP 鬆弛對比（流量變數連續、W/S 變 0..1 連續）
python gurobi/sc_model_gurobi.py --relax

# 設定 60 秒時間上限與求解策略
python gurobi/sc_model_gurobi.py --timelimit 60 --focus 1
```

輸出會包含：最佳目標值、MIPGap、開啟之倉庫/門市、主要非零流量與生產量等，用來和 C++ 版結果做對照。

---

### 4) 小提醒（hpclogin）

* 某些學校登入節點不鼓勵重度計算；若有 Slurm/作業系統，請改用計算節點（例如 `srun`/`sbatch`）執行。
* 若 WLS 因防火牆無法對外連線，請改用 **Named-User** 授權（需要 `grbgetkey`；可用 `conda install -c gurobi gurobi` 安裝授權工具後在**同一台機器**簽發），或改在你本機/WSL 跑 Gurobi。

---

### 5) C++ 與 Gurobi 的對解比對建議

* **相同參數來源**：`sc_params.hpp` 與 `gurobi/sc_params.py` 已對齊（集合/維度/順序一致）。
* **一致的模型**：`sc_model.cpp` 與 `gurobi/sc_model_gurobi.py` 使用同一組目標與限制（含 Big-M 與容量的體積換算）。
* 先跑 `--relax` 檢查 **LP 鬆弛** 是否一致，再跑整數版對照 **最優值與解構形**（開倉/開店、主要流量）。

> 有需要我可以加一個小比對腳本，把 C++ 解與 Gurobi 解整理成同一張表（含 top-k 流量/產量差異、目標值差、gap）做報告用。
