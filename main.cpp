#include <iostream> // std::flush
#include <stdio.h> // printf
#include <cstdint> // uint32_t
#include <cmath>
#include <string>
#include <algorithm> // for std::sort
#include <limits> // fp64 inf, nan
#include <vector>
#include <unordered_map> // map
#include <queue> // min-heap
#include <chrono> // 測速
#include <optional> // for node queue parallel

using namespace std;

class FOP {
public:
	static constexpr double EPS = 1e-4; // 浮點誤差
	
	static bool isInt(double x) { // 浮點數 x 是否是整數
		return abs(x - round(x)) <= EPS;
	}
	
	static bool isZero(double x) { // 浮點數 x 是否為 0
		return abs(x) <= EPS;
	}
	
	static bool isPos(double x) { // 浮點數 x 是否為正數
		return x >= EPS;
	}
};

// 平行化區域 start
bool enableMatrixEliminationParallel = false; // 啟用矩陣列運算 avx256 向量化加速

// 平行化版本的 array elimination, 用 A_{ij} 消去行 j 的其他元素, 並將列 i 同除 A_{ij}, 使 A_{ij} = 1
#include <omp.h>
#include <immintrin.h>
void parallelArrayElimination(uint32_t cols, vector<double>& arr, uint32_t i, uint32_t j) { // arr 是扁平化的二維陣列
	uint32_t rows = arr.size() / cols;
	
	const uint32_t simdWidth = 4; // 一次處理 4 個 double
	double* ptrRowI = arr.data() + cols * (size_t)(i); // A_{i0} 的指標
	
	for (uint32_t k = 0; k < rows; k++) { // k 為每個 pthread 負責的列運算 row 編號, 將列 i 乘常數消去列 k
		if (k == i) continue; // 列 i 不能消去列 i
		
		double* ptrRowK = arr.data() + cols * (size_t)(k); // A_{k0} 的指標
		const double ratio = ptrRowK[j] / ptrRowI[j]; // 列 i 乘常數加到列 k, 消去行 j 的其他元素
		// if (FOP::isZero(ratio)) continue; // [我不知道為什麼加這一行會卡死] 常數是 0 就不做列運算
		
		__m256d vecRatio = _mm256_set1_pd(ratio); // 全部都是常數的向量
		uint32_t c = 0;
		for (; c + (simdWidth-1) < cols; c += simdWidth) { // 向量運算不能超出列 k, 不然會寫到列 k+1 的記憶體
			__m256d vecScaled = _mm256_mul_pd(_mm256_loadu_pd(ptrRowI + c), vecRatio); // 列 i 乘常數 (vec)
			_mm256_storeu_pd(ptrRowK + c, _mm256_sub_pd(_mm256_loadu_pd(ptrRowK + c), vecScaled)); // 列 i 乘常數消去列 k (vec)
		}
		for (; c < cols; c++) ptrRowK[c] -= ptrRowI[c] * ratio; // 列 i 乘常數消去列 k
		
		ptrRowK[j] = 0; // 保證行 j 的其他元素被消去
	}
	
	const double aij = arr[cols * i + j];
	for (uint32_t k = 0; k < cols; k++) arr[cols * i + k] /= aij; // 列 i 同除 A_{ij}, 使 A_{ij} = 1
}
// 平行化區域 end

// 以下為單執行序的原始演算法

const double FP64_INF = numeric_limits<double>::infinity();
const double FP64_NAN = numeric_limits<double>::quiet_NaN();

enum class Relation { LEQ, EQ, GEQ }; // <= (LEQ), = (EQ), >= (GEQ)

template <typename K, typename V>
bool mapHasKey(const unordered_map<K, V>& m, const K& key) {
	return m.find(key) != m.end(); // 如果 unordered_map 存在某個 key, 則回傳 value, 否則回傳 defaultValue
}

class VarBimap { // 變數名稱與編號的雙向映射
private:
	unordered_map<string, uint32_t> strToIndex;
	vector<string> indexToStr;

public:
	uint32_t getVarCount() const { // 獲取已註冊變數的數量
		return indexToStr.size();
	}
	
	uint32_t getVarIndex(const string& varName) { // 註冊一個字串變數, 自動分配編號, 回傳這個字串變數對應的編號 j (x_j)
		if (mapHasKey<string, uint32_t>(strToIndex, varName)) return strToIndex[varName]; // 如果字串變數已註冊過, 回傳對應的編號 j
		
		const uint32_t newVarIndex = getVarCount(); // 如果字串變數沒有註冊過, 分配編號 0, 1, 2, ...
		strToIndex[varName] = newVarIndex; // 註冊雙向映射
		indexToStr.push_back(varName);
		return newVarIndex; // 回傳分配的新編號
	}
	
	string getVarName(uint32_t varIndex) const { // 編號 j (x_j) 轉字串變數
		if (varIndex <= indexToStr.size() - 1) return indexToStr[varIndex];
		return "[unknown-var]";
	}
};

class Linearform { // 線性函數
public:
	unordered_map<uint32_t, double> terms; // x_j index 映射到係數
	
	Linearform& add(double coef, uint32_t varIndex) { // 添加一個未知數項 (chaining)
		if (mapHasKey<uint32_t, double>(terms, varIndex)) terms[varIndex] += coef; // x_j 存在就加上係數
		else terms[varIndex] = coef; // x_j 不存在就建立未知數項
		
		return *this;
	}
	
	void negate() { // 對這個線性函數 *-1
		for (auto& term: terms) term.second = -term.second; // 對所有係數變號
	}
	
	void print(VarBimap& bimap, bool addNewLine = true) { // debug
		uint32_t c = 0;
		for (uint32_t i = 0; i < bimap.getVarCount(); i++) if (mapHasKey<uint32_t, double>(terms, i)) {
			if (c++) printf(" + ");
			printf("%.2f[%s]", terms[i], bimap.getVarName(i).c_str());
		}
		if (addNewLine) printf("\n");
	}
};

class Constraint { // 約束. 由一個線性函數, 關係, 和右側常數組成: c1 x1 + c2 x2 + ... ~ r
private:
	Linearform linearform; // 線性函數
	Relation relation = Relation::EQ; // >= or = or <=
	double rightConst = 0; // 右側常數
	
	void set(Relation relation, double rightConst) { // leq, eq, geq 共用的 set 邏輯
		this->relation = relation;
		this->rightConst = rightConst;
	}

public:
	Constraint& add(double coef, uint32_t varIndex) { // 添加一個未知數項 (chaining)
		linearform.add(coef, varIndex);
		return *this;
	}
	
	Constraint& leq(double rightConst) { // 添加最右側的 "<= r" 部分 (chaining)
		set(Relation::LEQ, rightConst);
		return *this;
	}
	
	Constraint& eq(double rightConst) { // 添加最右側的 "= r" 部分 (chaining)
		set(Relation::EQ, rightConst);
		return *this;
	}
	
	Constraint& geq(double rightConst) { // 添加最右側的 ">= r" 部分 (chaining)
		set(Relation::GEQ, rightConst);
		return *this;
	}
	
	unordered_map<uint32_t, double> getTerms() {
		return linearform.terms;
	}
	
	double getRightConst() {
		return rightConst;
	}
	
	void stdOfNegativeRightConst() { // 對負的右側常數進行標準化
		if (rightConst >= 0) return; // 若右側常數 >= 0, 跳過這一步
		
		rightConst = -rightConst; // 對右側常數變號
		linearform.negate(); // 對所有係數變號
		
		if (relation == Relation::LEQ) relation = Relation::GEQ; // 兩側變號需要將 <= 和 >= 轉向
		else if (relation == Relation::GEQ) relation = Relation::LEQ;
	}
	
	bool haveSlackVar() { // 是否需要添加 slack var
		return relation != Relation::EQ; // 只有 >= 和 <= 會產生 slack var
	}
	
	double getSlackVarCoef() { // 如果有 slack var, 回傳係數, 否則回傳 0
		if (relation == Relation::LEQ) return 1; // "... <= r" -> "... + x = r"
		if (relation == Relation::GEQ) return -1; // "... >= r" -> "... - x = r"
		return 0; // "=" 沒有 slack var
	}
	
	bool haveArtificialVar() { // 是否需要添加 artificial var
		return relation != Relation::LEQ; // 只有 = 和 >= 會產生 artificial var
	}
	
	void print(VarBimap& bimap, bool addNewLine = true) { // debug
		linearform.print(bimap, false);
		
		if (relation == Relation::LEQ) printf(" <= ");
		else if (relation == Relation::EQ) printf(" = ");
		else if (relation == Relation::GEQ) printf(" >= ");
		
		printf("%.2f", rightConst);
		
		if (addNewLine) printf("\n");
	}
};

class LP { // Linear Programming
private:
	class Tableau { // tableau, 注意: 第零列為 head, 第一列開始才是約束的部分
	private:
		vector<double> arr; // 扁平化的二維陣列
		
		double& arr_(uint32_t i, uint32_t j) { // 訪問扁平化的二維陣列
			return arr[cols * i + j];
		}
	
	public:
		uint32_t rows; // row 數
		uint32_t cols; // col 數
		vector<int32_t> baseVarIndexs; // 基底變數的編號, -1 代表 artifical var
		
		double& operator()(uint32_t i, uint32_t j) { // 訪問扁平化的二維陣列
			return arr[cols * i + j];
		};
		
		void init(uint32_t rows, uint32_t cols) { // 初始化 tableau
			this->rows = rows; // 設定列數
			this->cols = cols; // 設定行數
			arr = vector<double>(rows * cols, 0); // 預設全為零
			baseVarIndexs = vector<int32_t>(rows, 0); // 因為第零列沒有基底編號, 使 index 對齊
		}
		
		void scaleRow(uint32_t i, double s) { // 將列 i 除以常數 s
			for (uint32_t k = 0; k < cols; k++) arr_(i, k) /= s;
		}
		
		void addRowToRow(uint32_t i, uint32_t j, double s) { // 列 i 乘常數 s 加到列 j
			for (uint32_t k = 0; k < cols; k++) arr_(j, k) += arr_(i, k) * s;
		}
		
		void elimination(uint32_t i, uint32_t j) { // 用 A_{ij} 消去行 j 的其他元素, 並將列 i 同除 A_{ij}, 使 A_{ij} = 1
			if (enableMatrixEliminationParallel) { // 啟用矩陣列運算 avx256 向量化加速
				parallelArrayElimination(cols, arr, i, j);
				return; // 跳過原始演算法
			}
			
			const double aij = arr_(i, j);
			for (uint32_t k = 0; k < rows; k++) {
				if (k != i) addRowToRow(i, k, -arr_(k, j) / aij); // 用 A_{ij} 消去行 j 的其他元素
			}
			// 消去完成後,再統一設為 0
			for (uint32_t k = 0; k < rows; k++) {
				if (k != i) arr_(k, j) = 0;
			}
			scaleRow(i, aij); // 將列 i 同除 A_{ij}, 使 A_{ij} = 1
		}
		
		void print() { // debug
			for (uint32_t i = 0; i < rows; i++) {
				for (uint32_t j = 0; j < cols; j++) printf("%.02f ", (*this)(i, j));
				if (i > 0) printf("= x_%d", baseVarIndexs[i]); // print 基底
				printf("\n");
			}
		}
	} tableau;
	
	bool isMin; // min/max
	Linearform& objFunc; // 目標函數
	vector<Constraint>& multiCon; // 多個約束
	
	uint32_t varCount; // 一般變數的個數, index 為 0 ~ varCount-1
	vector<Constraint> varRangeMultiCon; // 將變數範圍轉為多個約束
	
	int32_t findNewBaseVarIndex() { // 尋找一個新的基底變數, 若沒找到則回傳 -1
		for (uint32_t j = 0; j <= tableau.cols - 2; j++) if (FOP::isPos(tableau(0, j))) return j; // 最後一列是基底常數, 不能進入
		return -1;
	}
	
	int32_t findMinPosRatioRowIndex(uint32_t baseVarIndex) { // 選定要進入的基底後, 尋找一個 Aij / r 最小的正比值, 回傳這個值在第幾列, 找不到回傳 -1
		double minPosRatio = 1e300; // 最小正比值: 右側常數/係數
		int32_t minPosRatioRowIndex = -1; // 要更換基底的 row index, 若為 -1 代表沒有找到
		for (uint32_t i = 1; i < tableau.rows; i++) if (FOP::isPos(tableau(i, baseVarIndex))) { // 要正比值
			double ratio = tableau(i, tableau.cols - 1) / tableau(i, baseVarIndex);
			if (ratio < minPosRatio) {
				minPosRatio = ratio;
				minPosRatioRowIndex = i;
			}
		}
		return minPosRatioRowIndex;
	}
	
	bool isTableauHaveArtificialVar() { // tableau 的基底是否含有 artificial var
		for (uint32_t i = 1; i < tableau.rows; i++) if (tableau.baseVarIndexs[i] == -1) return true;
		return false;
	}
	
	void handleUnbound(uint32_t newBaseVarIndex) { // 若無界, 輸出當前向量和一個無限增長的方向向量 (這個向量集合仍然符合所有約束)
		solutionType = Type::UNBOUNDED; // 無界
		
		solution = vector<double>(varCount, 0); // 空的解向量, 長度為一般變數的數目
		unboundedDirection = vector<double>(varCount, 0); // 方向向量
		for (uint32_t i = 1; i < tableau.rows; i++) {
			uint32_t baseVarIndex = tableau.baseVarIndexs[i];
			if (baseVarIndex <= varCount - 1) { // 如果基底變數編號為 0 ~ varCount-1 代表為一般變數
				solution[baseVarIndex] = tableau(i, tableau.cols - 1);
				unboundedDirection[baseVarIndex] = tableau(i, newBaseVarIndex) * (isMin ? 1 : -1);
			}
		} // slack var 編號 >= varCount, 所以不會出現在解向量裡
		
		extremum = isMin ? -FP64_INF : FP64_INF; // min -> -inf ; max -> inf
	}
	
	bool runMinSimplexMethod() { // 對 tableau 執行 min simplex method, 若有界回傳 true, 無解或無界回傳 false
		while (true) {
			const int32_t newBaseVarIndex = findNewBaseVarIndex(); // 嘗試尋找新基底 [複雜度: n]
			if (newBaseVarIndex == -1) break; // 若沒有找到可進入的基底, 跳出迴圈
			
			const int32_t rowIndex = findMinPosRatioRowIndex(newBaseVarIndex); // 嘗試尋找最小正數比值的列編號 [複雜度: m]
			if (rowIndex == -1) { // 若新基底存在, 但最小正數比值不存在, 則 LP 問題無界
				handleUnbound(newBaseVarIndex);
				return false; // 提前結束迴圈
			}
			
			tableau.elimination(rowIndex, newBaseVarIndex); // 對新基底的行做消元, 只留下最小正數比值的列 [複雜度: m*n]
			tableau.baseVarIndexs[rowIndex] = newBaseVarIndex; // 更改新基底編號
		}
		
		if (isTableauHaveArtificialVar()) return false; // 執行完 simplex method, 若最小的 L1-norm 起始向量包含 artificial var, 則 LP 問題無解
		return true;
	}
	
	void initTableau() { // init tableau
		uint32_t slackVarCount = 0; // 計算 slack var 個數, 決定 tableau 的 col 數 (因為要分配連續記憶體)
		for (Constraint& con: multiCon) if (con.haveSlackVar()) slackVarCount++;
		for (Constraint& con: varRangeMultiCon) if (con.haveSlackVar()) slackVarCount++;
		tableau.init(1 + multiCon.size() + varRangeMultiCon.size(), varCount + slackVarCount + 1); // init tableau
	}
	
	void insertConToTableau() { // 將約束插入 tableau
		uint32_t rowIndex = 1;
		uint32_t slackVarColIndex = varCount - 1; // 因為一般變數的 col index 為 0 ~ varCount-1, 所以 slack var 插入的 col index 從這裡開始數
		auto setTableauRow = [&](Constraint& con) { // 將一個約束加入到 tableau
			for (auto& [varIndex, coef]: con.getTerms()) tableau(rowIndex, varIndex) = coef; // 填入一般變數
			
			if (con.haveSlackVar()) tableau(rowIndex, ++slackVarColIndex) = con.getSlackVarCoef(); // 如果有 slack var, 需要添加係數 1 或 -1 到 tableau 內
			
			tableau(rowIndex, tableau.cols - 1) = con.getRightConst(); // 在列的最右元素, 設定右側常數 (基底的值)
			tableau.baseVarIndexs[rowIndex] = con.haveArtificialVar() ? -1 : slackVarColIndex; // -1 代表 artifical var
			
			rowIndex++;
		};
		for (Constraint& con: multiCon) setTableauRow(con);
		for (Constraint& con: varRangeMultiCon) setTableauRow(con);
	}
	
	bool checkInfeasible() { // 檢查是否無解, 並做處理, 若無解則回傳 false, 非無解則回傳 true
		if (!isTableauHaveArtificialVar()) return true; // 如果不存在 artificial var, 跳過這一步
		
		for (uint32_t i = 1; i < tableau.rows; i++) if (tableau.baseVarIndexs[i] == -1) {
			tableau.addRowToRow(i, 0, 1); // 將 artificial var 的列加到第零列, 消去第零列的 art-var 係數 -1 (但我們的演算法不會儲存 art-var 的行係數)
		}
		
		if (!runMinSimplexMethod()) return false; // 嘗試求出一個最小的 L1-norm 起始向量, 若無解則回傳 false
		
		for (uint32_t j = 0; j < tableau.cols; j++) tableau(0, j) = 0; // 雖然理論上第零列應該是全 0 的, 只是為了消除小誤差
		return true; // 可能為有界或無界
	}
	
	void insertObjFuncToTableau() { // 將目標函數插入到 tableau 的第零列
		for (auto& [varIndex, coef]: objFunc.terms) tableau(0, varIndex) = coef * (isMin ? -1 : 1);
	} // 此時 tableau 的第零列是空的, 填入目標函數. 注意: 因為將 max obj func 變號會轉為 min 問題, 所以 max 問題在這裡會填入 +coef 而不是 -coef
	
	void handleNonZeroHead() { // 處理 row 0 的非零基底行元素 (phase-2). 若沒做 phase-1 則這一步不會有實際影響
		for (uint32_t i = 1; i < tableau.rows; i++) {
			uint32_t baseVarIndex = tableau.baseVarIndexs[i]; // 檢查每一列對應的基底變數所對應的行的 row 0 元素是不是 0
			if (!FOP::isZero(tableau(0, baseVarIndex))) tableau.addRowToRow(i, 0, -tableau(0, baseVarIndex)); // 如果非 0, 需要執行列運算, 消去它
		}
	}
	
	void handleBound() { // 處理有界的情況
		solutionType = Type::BOUNDED;
		
		solution = vector<double>(varCount, 0); // 空的解向量, 長度為一般變數的數目
		for (uint32_t i = 1; i < tableau.rows; i++) {
			uint32_t baseVarIndex = tableau.baseVarIndexs[i];
			if (baseVarIndex <= varCount - 1) solution[baseVarIndex] = tableau(i, tableau.cols - 1); // 如果基底變數編號為 0 ~ varCount-1 代表為一般變數
		} // slack var 編號 >= varCount, 所以不會出現在解向量裡
		
		extremum = tableau(0, tableau.cols - 1) * (isMin ? 1 : -1); // 極值
	}

public:
	enum class Type { BOUNDED, UNBOUNDED, INFEASIBLE }; // 有界, 無界, 無解
	
	Type solutionType; // 解的狀態
	vector<double> solution; // 解向量, 若無解會是空的
	vector<double> unboundedDirection; // 若無界, 此值會是一個方向向量, 即使往無窮遠移動仍然滿足目標函數
	double extremum; // min/max 極值
	
	LP(bool isMin, Linearform& objFunc, vector<Constraint>& multiCon, vector<pair<double, double>>& varRange)
	: isMin(isMin), objFunc(objFunc), multiCon(multiCon) {
		varCount = varRange.size();
		uint32_t i = 0;
		for (auto& [varMin, varMax]: varRange) { // min <= x_i <= max, 將變數範圍轉為約束
			if (varMin > 0) varRangeMultiCon.push_back(Constraint().add(1, i).geq(varMin)); // x_i >= min
			if (!isinf(varMax)) varRangeMultiCon.push_back(Constraint().add(1, i).leq(varMax)); // x_i <= max
			i++;
		}
		
		initTableau(); // init tableau
		insertConToTableau(); // 將約束插入 tableau
		
		if (!checkInfeasible()) { // 檢查是否無解, 並做處理 (phase-1)
			solutionType = Type::INFEASIBLE; // checkInfeasible 輸出 false, 無解
			extremum = FP64_NAN;
		} else { // checkInfeasible 輸出 true, 可能為有界或無界. checkInfeasible 已處理基底不足的問題
			insertObjFuncToTableau(); // 將目標函數插入到 tableau 的第零列
			handleNonZeroHead(); // 處理 row 0 的非零基底行元素 (phase-2)
			if (runMinSimplexMethod()) handleBound(); // 處理有界情況, 如果無界會提前終止並處理
		}
	}
	
	void print(bool showCon = false) { // debug
		if (showCon) {
			VarBimap bimap; // 因為 IP 已經將字串變數轉為 index 跟 LP 溝通, 所以 LP 抽象層並不知道 varName, 所以這邊註冊一個 x0, x1, x2 (抽象 index)
			for (uint32_t i = 0; i < varCount; i++) bimap.getVarIndex("x" + to_string(i)); // 註冊 x0, x1, x2, ...
			
			printf(isMin ? "min " : "max "); // 印出目標函數
			objFunc.print(bimap);
			for (Constraint& con: multiCon) con.print(bimap); // 印出約束
			
			if (varRangeMultiCon.size() == 0) printf("Var range is empty.");
			for (Constraint& con: varRangeMultiCon) { // 印出變數範圍
				con.print(bimap, false);
				printf("; ");
			}
			printf("\n");
		}
		
		printf("Type: "); // 印出 LP 的解的類型
		if (solutionType == Type::BOUNDED) printf("Bounded\n");
		else if (solutionType == Type::UNBOUNDED) printf("Unbounded\n");
		else if (solutionType == Type::INFEASIBLE) printf("Infeasible\n");
		
		printf(isMin ? "LP Minimum = " : "LP Maximum = "); // 印出極值
		printf("%.2f\n", extremum);
		
		printf("LP Solution: "); // 印出解向量
		for (size_t i = 0; i < solution.size(); i++) printf("x%d = %.2f; ", (int)i, solution[i]);
		printf("\n");
		
		if (unboundedDirection.size() > 0) { // 如果無界, 印出方向向量
			printf("Unbounded delta: ");
			for (size_t i = 0; i < unboundedDirection.size(); i++) printf("x%d = %.2f; ", (int)i, unboundedDirection[i]);
			printf("\n");
		}
	}
};

class IP { // Integer Programming
private:
	class Node { // branch & bound 的 node
	private:
		static int32_t getSplitVarIndex(LP& lp) { // 決定下次分支要切分的基底變數編號, 如果基底變數的值全部都是整數 (LP feasible) 則回傳 -1
			for (uint32_t i = 0; i < lp.solution.size(); i++) if (!FOP::isInt(lp.solution[i])) return i; // 目前的策略是尋找編號最小的 float 變數值
			return -1;
		}
		
		static double getSplitValue(pair<double, double> varRange, double varSolution) { // 根據某個基底變數的範圍和 LP 解, 決定切分值
			return varSolution; // 後來發現直接這樣切會更快
			
			if (varRange.second == FP64_INF) return varSolution; // 將 [?, inf] 切分為 [?, sol] & [sol, inf]
			return (varRange.first + varRange.second) / 2; // 將 [a, b] 根據中點 (a+b)/2 切分
		} // 因為如果照著某個 LP solution 去分支 (投影片的演算法), 會嘗試對特定變數切 1000, 999, 998, ... (所以目前是直接對半切)
	
	public:
		struct cmp { // priority queue 用的比較子, 會變成以 node 下界排序的 min-heap
			bool operator()(const Node& a, const Node& b) const { return a.lowerBound > b.lowerBound; }
		};
		
		enum class Type { IP_FEASIBLE, LP_FEASIBLE, INFEASIBLE, UNBOUNDED }; // 有 IP 解, 有 LP 解, 無解, 無界 (比較麻煩)
		
		vector<double> solution; // float LP 解
		double lowerBound = FP64_NAN; // 節點的 float min LP 下界, int 解只可能更大
		Type type; // node 的型態
		vector<pair<double, double>> varRangeLeft; // 下一次分支, 左子節點每個變數的範圍, 注意: 變數 x_j 的範圍為 [varRange[j].first, varRange[j].second]
		vector<pair<double, double>> varRangeRight; // 下一次分支, 右子節點每個變數的範圍
		
		Node(Linearform& objFunc, vector<Constraint>& multiCon, vector<pair<double, double>>& varRange) {
			LP lp(true, objFunc, multiCon, varRange); // 解 LP (已經將 max 標準化為 min)
			solution = lp.solution; // float LP 解
			lowerBound = lp.extremum; // float min LP 的極值
			
			if (lp.solutionType == LP::Type::INFEASIBLE) type = Type::INFEASIBLE;
			else if (lp.solutionType == LP::Type::UNBOUNDED) type = Type::UNBOUNDED;
			else if (lp.solutionType == LP::Type::BOUNDED) {
				int32_t splitVarIndex = getSplitVarIndex(lp); // 下次分支要切分的基底變數編號
				if (splitVarIndex == -1) { // 如果基底變數的值全部都是整數 (IP feasible), 更新全域的 int min IP 上界
					type = Type::IP_FEASIBLE;
				} else { // 如果有基底變數的值是 float, 計算出下一次分支的左右節點的變數範圍
					type = Type::LP_FEASIBLE;
					varRangeLeft = varRange; // 複製兩份 varRange, 分別為左右子節點的變數範圍
					varRangeRight = varRange;
					
					double splitValue = floor(getSplitValue(varRangeLeft[splitVarIndex], lp.solution[splitVarIndex])); // 切分值
					varRangeLeft[splitVarIndex].second = splitValue; // 將左子節點的切分基底值的上界設為 splitValue
					varRangeRight[splitVarIndex].first = splitValue + 1; // 將右子節點的切分基底值的下界設為 splitValue + 1
				}
			}
		}
		
		void print(VarBimap bimap, bool showRange = false) {
			if (type == Type::IP_FEASIBLE) {
				printf("[IP solution found] Objective value = %.2f\n", lowerBound);
			} else if (type == Type::INFEASIBLE) {
				printf("IP solution is infeasible.\n");
			} else if (type == Type::UNBOUNDED) {
				printf("LP solution is unbounded, maybe IP solution does not exist.\n");
			} else if (type == Type::LP_FEASIBLE) {
				printf("IP feasible objective value >= %.2f\n", lowerBound);

				// 只有 LP_FEASIBLE 才會有左右子區間；先檢查向量是否為空
				if (!varRangeLeft.empty() && !varRangeRight.empty()) {
					for (size_t i = 0; i < solution.size(); ++i) {
						if (varRangeLeft[i].first != varRangeRight[i].first) {
							printf("Next split: %s = ", bimap.getVarName((uint32_t)i).c_str());
							printf("[%.2f %.2f] & ", varRangeLeft[i].first, varRangeLeft[i].second);
							printf("[%.2f %.2f]\n", varRangeRight[i].first, varRangeRight[i].second);
							break;
						}
					}
				}
			}

			if (showRange && type == Type::LP_FEASIBLE && !varRangeLeft.empty() && !varRangeRight.empty()) {
				printf("Left child node var range: ");
				for (size_t vi = 0; vi < varRangeLeft.size(); ++vi) {
					const auto &rng = varRangeLeft[vi];
					if (rng.first != 0 || rng.second != FP64_INF) {
						printf("%s = [%.2f, %.2f]; ", bimap.getVarName((uint32_t)vi).c_str(), rng.first, rng.second);
					}
				}
				printf("\n");

				printf("Right child node var range: ");
				for (size_t vi = 0; vi < varRangeRight.size(); ++vi) {
					const auto &rng = varRangeRight[vi];
					if (rng.first != 0 || rng.second != FP64_INF) {
						printf("%s = [%.2f, %.2f]; ", bimap.getVarName((uint32_t)vi).c_str(), rng.first, rng.second);
					}
				}
				printf("\n");

			}
		}
	};
	
	bool isMin; // min = 1, max = 0
	Linearform objFunc; // 目標函數
	vector<Constraint> multiCon; // 多個約束
	
	VarBimap bimap; // 變數映射
	priority_queue<Node, vector<Node>, Node::cmp> nodeQueue; // 以 float LP 下界排序的 min-heap, 先展開下界較小的 node 比較容易找到更小的解
	double objValueUpperBound = FP64_INF; // 因為是求 min IP 問題, 所以有一個全域上界
	
	uint32_t nodeSolvedCount = 0; // [debug 變數] 計算了幾次 LP 問題
	int64_t lastTimePrintNodeInfo = getSystemTimeSec(); // [debug 變數] 上一次印出 node queue 資訊的時間
	
	void init() { // 初始化 IP 問題
		if (!isMin) objFunc.negate(); // 將 max 問題轉為 min 問題, 只需要將目標函數變號即可
		
		uint32_t varCount = bimap.getVarCount(); // 一般變數的個數
		vector<pair<double, double>> varRange(varCount, { 0, FP64_INF }); // 生成一般變數的範圍, branch & bound 的 root node 的變數範圍全為 [0, inf]
		Node rootNode = Node(objFunc, multiCon, varRange); // root node
		checkNode(rootNode); // 檢查 node 的 solution type
	}
	
	void checkNode(Node& node) { // 檢查一個 node 的 solution type, 決定是否要更新全域上界或推入 min heap
		if (node.type == Node::Type::IP_FEASIBLE && node.lowerBound < objValueUpperBound) {
			solutionType = Type::BOUNDED;
			solution = node.solution; // 更新全域 IP 解
			objValueUpperBound = node.lowerBound; // [剪枝] 如果 node 有整數解向量, 並且比現有的解更好, 更新全域上界, 不用繼續往下尋找
		} // [剪枝] 如果 node 有整數解向量, 沒有比現有的解更好, 無視
		else if (node.type == Node::Type::LP_FEASIBLE && node.lowerBound < objValueUpperBound) { // 如果 node 有浮點解向量
			nodeQueue.push(node); // 將 root node push 進 min-heap (繼續往下搜尋)
		} // [剪枝] "node 下界 >= 全域上界" 的分支不用繼續往下搜尋, 因為無法取得更好的結果
		else if (node.type == Node::Type::UNBOUNDED) { // 如果 node 無界
			solutionType = Type::UNBOUNDED; // 停止計算 IP
		} // [剪枝] 如果 node 無解, 無視它
		
		nodeSolvedCount++; // 解 LP 式子的次數與 check 次數相同
		
		return; // [testing] 目前禁用 print node info
		const int64_t systemTimeNowSec = getSystemTimeSec();
		if (systemTimeNowSec != lastTimePrintNodeInfo) { // 如果時間戳改變, print 一次 node 資訊
			lastTimePrintNodeInfo = systemTimeNowSec;
			printf("[ Node queue size = %d | %d LP nodes solved ]\n", (uint32_t)nodeQueue.size(), nodeSolvedCount);
		}
	}
	
	int64_t getSystemTimeSec() { // 獲取目前系統時間戳 (sec)
		const auto epochTime = chrono::system_clock::now().time_since_epoch();
		return chrono::duration_cast<chrono::seconds>(epochTime).count();
	}

public:
	enum class Type { BOUNDED, INFEASIBLE, UNBOUNDED };
	
	Type solutionType = Type::INFEASIBLE; // 預設是無解, 如果有發現 IP 解會修改此值
	vector<double> solution; // 全域 IP 解向量
	double extremum; // min/max 極值
	
	IP(const string& mode, vector<pair<double, string>> terms) { // 宣告 min/max 和目標函數
		isMin = mode == "min"; // min/max
		for (auto& [coef, varName]: terms) objFunc.add(coef, bimap.getVarIndex(varName)); // 目標函數
	}
	
	IP& addConstraint(vector<pair<double, string>> terms, const string& mode, double rightConst) { // 添加約束 (chaining)
		Constraint con;
		for (auto& [coef, varName]: terms) con.add(coef, bimap.getVarIndex(varName)); // 添加目標函數
		
		if (mode == "<=") con.leq(rightConst); // 添加右側常數
		else if (mode == ">=") con.geq(rightConst);
		else con.eq(rightConst);
		
		con.stdOfNegativeRightConst(); // 對負的右側常數進行標準化
		multiCon.push_back(con);
		
		return *this;
	}
	
	void solve() { // 計算 IP 問題
		init(); // 生成初始 node 並 push 進 min-heap
		
		while (nodeQueue.size() > 0 && solutionType != Type::UNBOUNDED) { // min-heap 還有 node 就繼續分支, 目前 unbounded 會強制停下
			Node node = nodeQueue.top(); // 取出下界較小的 node 比較容易找到更小的解
			nodeQueue.pop();
			
			Node leftChildNode(objFunc, multiCon, node.varRangeLeft); // 生成並計算左子節點的 LP 問題
			Node rightChildNode(objFunc, multiCon, node.varRangeRight); // 生成並計算右子節點的 LP 問題
			checkNode(leftChildNode); // 檢查 child node 的解
			checkNode(rightChildNode);
		}
		
		extremum = objValueUpperBound * (isMin ? 1 : -1); // 極值
	}
	
	void solveParallel() { // 計算 IP 問題 (node level parallel)
		init(); // 生成初始 node 並 push 進 min-heap
		
		uint32_t workingThreadCount = 0; // 正在計算中的 thread
		
		#pragma omp parallel // 建立一個執行緒池
		{
			optional<Node> nodeOpt;
			
			while (true) { // 不斷 busy waiting
				if (solutionType == Type::UNBOUNDED) break; // 如果有 node 的 LP 解出現 unbounded 會強制停下
				
				bool hasNode = false;
				bool hasActiveThread = false;
				#pragma omp critical // 同一時間只有一個執行緒可以嘗試從 queue 取出 node
				{
					while (nodeQueue.size() > 0) {
						nodeOpt = nodeQueue.top(); // 取出下界較小的 node 比較容易找到更小的解
						nodeQueue.pop();
						
						if (nodeOpt.value().lowerBound < objValueUpperBound) { // 在 critical section 不斷嘗試取出一個 "node 下界 < 全域上界" 的 node
							workingThreadCount++;
							hasNode = true;
							break;
						}
					}
					if (workingThreadCount > 0) hasActiveThread = true;
				}
				if (!hasNode && hasActiveThread) continue; // node queue 是空的, 但還有 thread 還在計算, 就 busy waiting
				if (!hasNode && !hasActiveThread) break; // node queue 是空的, 並且沒有 thread 在計算, 代表 node tree 已遍歷完畢
				
				// 每個執行緒獨立計算自己的 LP 子問題
				Node& node = nodeOpt.value();
				Node leftChildNode(objFunc, multiCon, node.varRangeLeft); // 計算左子樹 (left node)
				Node rightChildNode(objFunc, multiCon, node.varRangeRight); // 計算右子樹 (right node)
				
				#pragma omp critical // 同一時間只有一個執行緒可以將 node 推入 queue, 修改 workingThreadCount 和 objValueUpperBound
				{
					checkNode(leftChildNode); // 檢查左右子樹的 node 的 LP 解, 決定是否要更新全域上界或推入 node queue
					checkNode(rightChildNode);
					workingThreadCount--; // 工作做完, 計數 -1
				}
			}
		}
		
		extremum = objValueUpperBound * (isMin ? 1 : -1); // 因為有將 max 問題轉為 min 問題, 極值要記得變號
	}
	
	uint32_t getNodeSolvedCount() {
		return nodeSolvedCount;
	}
	
	void print_grouped_solution(bool show_zero = false) const {
		struct Item { std::string name; double val; };
		std::vector<Item> items;
		items.reserve(solution.size());
		for (size_t i = 0; i < solution.size(); ++i) {
			items.push_back({ bimap.getVarName((uint32_t)i), solution[i] });
		}

		auto is_zero = [](double x){ return std::abs(x) <= FOP::EPS; };
		auto print_group = [&](const std::string& prefix, const char* title) {
			std::vector<Item> g;
			for (const auto& it : items) {
				// 以字首判斷群組：P[, X[, Y[, U[, W[, S[
				if (it.name.rfind(prefix, 0) == 0) {
					if (show_zero || !is_zero(it.val)) g.push_back(it);
				}
			}
			if (g.empty()) return;
			std::sort(g.begin(), g.end(), [](const Item& a, const Item& b){
				return a.name < b.name;
			});
			printf("%s\n", title);
			for (const auto& it : g) {
				if (FOP::isInt(it.val))
					printf("  %-28s = %.0f\n", it.name.c_str(), std::round(it.val));
				else
					printf("  %-28s = %.4f\n", it.name.c_str(), it.val);
			}
		};

		// 你可依需求調整順序
		print_group("W[", "Warehouses open (W_k):");
		print_group("S[", "Stores open (S_l):");
		print_group("P[", "Production (P[i,j]):");
		print_group("X[", "Shipments Factory→WH (X[i,j,k]):");
		print_group("Y[", "Shipments WH→Store (Y[i,k,l]):");
		print_group("U[", "Unmet demand (U[i,l]):");
	}
	
	void print(bool showCon = false, bool showGroupSolution = false) {
		if (showCon) {
			printf(isMin ? "min " : "max -> min ");
			objFunc.print(bimap);
			for (Constraint& con: multiCon) con.print(bimap);
		}

		printf("Type: ");
		if (solutionType == Type::BOUNDED) printf("Bounded\n");
		else if (solutionType == Type::UNBOUNDED) printf("Unbounded\n");
		else printf("Infeasible\n");

		printf(isMin ? "IP Minimum = " : "IP Maximum = ");
		printf("%.2f\n", extremum);

		printf("IP solution: ");
		for (size_t i = 0; i < solution.size(); i++) {
			printf("%s = %d; ", bimap.getVarName(i).c_str(), (int32_t)solution[i]);
		}
		printf("\n");

		printf("Number of LP nodes solved: %u\n", nodeSolvedCount);
		
		if (showGroupSolution) print_grouped_solution(true);
	}
};

#include "sc_params.hpp"
#include "sc_model.cpp"
class Tester { // 測速
private:
	int i, j, k, l;

public:
	Tester(int i, int j, int k, int l): i(i), j(j), k(k), l(l) {}
	
	pair<double, uint32_t> testOneIP(bool avx2MRO, bool nodeOmp) { // 測試單個 IP 問題的耗時和 LP node 解決數
		enableMatrixEliminationParallel = avx2MRO; // 啟用矩陣列運算 avx256 向量化加速
		bool enableNodeLevelParallel = nodeOmp; // node queue 會一次 pop 多個 node 做平行化計算
		
		SCParams P = default_sc_params(i, j, k, l); // 取參數 (可在 sc_params.hpp 改 default_sc_params() 內容)
		IP ip = build_supply_chain_ip(P); // 用參數建 IP 模型 (目標 + 限制)
		
		auto start = chrono::high_resolution_clock::now(); // 測速
		enableNodeLevelParallel ? ip.solveParallel() : ip.solve();
		auto end = chrono::high_resolution_clock::now();
		
		double exeTimeMs = chrono::duration<double, milli>(end - start).count(); // 執行總時間 (ms)
		return { exeTimeMs, ip.getNodeSolvedCount() };
	}
	
	pair<double, double> testParallel(uint32_t n, bool avx2MRO, bool nodeOmp) { // 測試 n 次同一個 IP 問題, 回傳平均耗時和平均 node 數
		printf("Solved IP problem count (avx2=%d omp=%d): ", avx2MRO, nodeOmp);
		cout << flush; // 分隔用
		
		double exeTimeMsSum = 0;
		uint32_t nodeSolvedCountSum = 0;
		for (uint32_t a = 0; a < n; a++) {
			auto [exeTimeMs, nodeSolvedCount] = testOneIP(avx2MRO, nodeOmp);
			exeTimeMsSum += exeTimeMs;
			nodeSolvedCountSum += nodeSolvedCount;
			cout << "*" << flush;
		}
		cout << endl;
		return { exeTimeMsSum / n, (double)nodeSolvedCountSum / n };
	}
	
	void test(uint32_t n) {
		auto [avgExeTimeMs_00, avgNodeSolvedCount] = testParallel(n, false, false);
		auto [avgExeTimeMs_10, _] = testParallel(n, true, false);
		auto [avgExeTimeMs_11, __] = testParallel(n, true, true);
		double avx2SpeedUp = avgExeTimeMs_00 / avgExeTimeMs_10; // avx2 matrix row operation speedup
		double ompSpeedUp = avgExeTimeMs_10 / avgExeTimeMs_11; // omp node level parallel speedup
		
		printf("-------------------- Tester --------------------\n");
		printf(" IP problem - Model parameters: (%d, %d, %d, %d)\n", i, j, k, l);
		printf(" Running %d IP problems\n", n);
		printf(" OpenMP max threads: %d\n", omp_get_max_threads());
		printf("------------------------------------------------\n");
		printf(" Average LP nodes solved per IP problem: %.0f\n", avgNodeSolvedCount);
		printf(" [AVX2: OFF, OMP: OFF] %.3f ms/IPprob\n", avgExeTimeMs_00);
		printf(
			" [AVX2: ON , OMP: OFF] %.3f ms/IPprob | AVX2 matrix row operation speedup: x %.2f (%.2f %%)\n",
			avgExeTimeMs_10, avx2SpeedUp, avx2SpeedUp / 4 * 100
		);
		printf(
			" [AVX2: ON , OMP: ON ] %.3f ms/IPprob | OpenMP node-level parallel speedup: x %.2f (%.2f %%)\n",
			avgExeTimeMs_11, ompSpeedUp, ompSpeedUp / omp_get_max_threads() * 100
		);
		printf("-------------------- Tester --------------------\n");
	}
};

int32_t main(int argc, char* argv[]) {
	Tester tester(3, 3, 3, 3);
	tester.test(10);
	
	return 0;
}
