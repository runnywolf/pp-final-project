#include <stdio.h> // printf
#include <cstdint> // uint32_t
#include <cmath>
#include <limits> // fp64 inf, nan
#include <vector>
#include <unordered_map> // map
#include <chrono> // 測速

using namespace std;

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
	uint32_t getVarCount() { // 獲取已註冊變數的數量
		return indexToStr.size();
	}
	
	uint32_t getVarIndex(const string& varName) { // 註冊一個字串變數, 自動分配編號, 回傳這個字串變數對應的編號 j (x_j)
		if (mapHasKey<string, uint32_t>(strToIndex, varName)) return strToIndex[varName]; // 如果字串變數已註冊過, 回傳對應的編號 j
		
		const uint32_t newVarIndex = getVarCount(); // 如果字串變數沒有註冊過, 分配編號 0, 1, 2, ...
		strToIndex[varName] = newVarIndex; // 註冊雙向映射
		indexToStr.push_back(varName);
		return newVarIndex; // 回傳分配的新編號
	}
	
	string getVarName(uint32_t varIndex) { // 編號 j (x_j) 轉字串變數
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
	static constexpr double EPS = 1e-10; // 浮點誤差
	
	enum class Type { BOUNDED, UNBOUNDED, INFEASIBLE }; // 有界, 無界, 無解
	
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
			baseVarIndexs = vector<int32_t>(cols, 0); // 因為第零列沒有基底編號, 使 index 對齊
		}
		
		void scaleRow(uint32_t i, double s) { // 將列 i 除以常數 s
			for (uint32_t k = 0; k < cols; k++) arr_(i, k) /= s;
		}
		
		void addRowToRow(uint32_t i, uint32_t j, double s) { // 列 i 乘常數 s 加到列 j
			for (uint32_t k = 0; k < cols; k++) arr_(j, k) += arr_(i, k) * s;
		}
		
		void elimination(uint32_t i, uint32_t j) { // 用 A_{ij} 消去行 j 的其他元素, 並將列 i 同除 A_{ij}, 使 A_{ij} = 1
			const double aij = arr_(i, j);
			for (uint32_t k = 0; k < rows; k++) if (k != i) addRowToRow(i, k, -arr_(k, j) / aij); // 用 A_{ij} 消去行 j 的其他元素
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
	
	Type solutionType; // 解的狀態
	vector<double> solution; // 解向量, 若無解會是空的
	vector<double> unboundedDirection; // 若無界, 此值會是一個方向向量, 即使往無窮遠移動仍然滿足目標函數
	double extremum; // min/max 極值
	
	int32_t findNewBaseVarIndex() { // 尋找一個新的基底變數, 若沒找到則回傳 -1
		for (uint32_t j = 0; j <= tableau.cols - 2; j++) if (tableau(0, j) >= EPS) return j; // 最後一列是基底常數, 不能進入
		return -1;
	}
	
	int32_t findMinPosRatioRowIndex(uint32_t baseVarIndex) { // 選定要進入的基底後, 尋找一個 Aij / r 最小的正比值, 回傳這個值在第幾列, 找不到回傳 -1
		double minPosRatio = 1e300; // 最小正比值: 右側常數/係數
		int32_t minPosRatioRowIndex = -1; // 要更換基底的 row index, 若為 -1 代表沒有找到
		for (uint32_t i = 1; i < tableau.rows; i++) if (tableau(i, baseVarIndex) >= EPS) { // 要正比值
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
			const int32_t newBaseVarIndex = findNewBaseVarIndex(); // 嘗試尋找新基底
			if (newBaseVarIndex == -1) break; // 若沒有找到可進入的基底, 跳出迴圈
			
			const int32_t rowIndex = findMinPosRatioRowIndex(newBaseVarIndex); // 嘗試尋找最小正數比值的列編號
			if (rowIndex == -1) { // 若新基底存在, 但最小正數比值不存在, 則 LP 問題無界
				handleUnbound(newBaseVarIndex);
				return false; // 提前結束迴圈
			}
			
			tableau.elimination(rowIndex, newBaseVarIndex); // 對新基底的行做消元, 只留下最小正數比值的列
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
		auto setTableauRow = [&](Constraint con) { // 將一個約束加入到 tableau
			for (auto& [varIndex, coef]: con.getTerms()) tableau(rowIndex, varIndex) = coef; // 填入一般變數
			
			if (con.haveSlackVar()) tableau(rowIndex, ++slackVarColIndex) = con.getSlackVarCoef(); // 如果有 slack var, 需要添加係數 1 或 -1 到 tableau 內
			
			tableau(rowIndex, tableau.cols - 1) = con.getRightConst(); // 在列的最右元素, 設定右側常數 (基底的值)
			tableau.baseVarIndexs[rowIndex] = con.haveArtificialVar() ? -1 : slackVarColIndex; // -1 代表 artifical var
			
			rowIndex++;
		};
		for (Constraint& con: multiCon) setTableauRow(con);
		for (Constraint& con: varRangeMultiCon) setTableauRow(con);
	}
	
	bool handleInfeasible() { // 檢查是否無解, 並做處理, 若無解則回傳 false, 非無解則回傳 true
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
			if (abs(tableau(0, baseVarIndex)) >= EPS) tableau.addRowToRow(i, 0, -tableau(0, baseVarIndex)); // 如果非 0, 需要執行列運算, 消去它
		} // 避免浮點誤差, 介於 +-EPS 內會視為 0
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
	LP(bool isMin, Linearform& objFunc, vector<Constraint>& multiCon, vector<pair<double, double>> varRange)
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
		
		if (!handleInfeasible()) { // 檢查是否無解, 並做處理 (phase-1)
			solutionType = Type::INFEASIBLE; // handleInfeasible 輸出 false, 無解
			extremum = FP64_NAN;
		} else { // handleInfeasible 輸出 true, 可能為有界或無界. handleInfeasible 已處理基底不足的問題
			insertObjFuncToTableau(); // 將目標函數插入到 tableau 的第零列
			handleNonZeroHead(); // 處理 row 0 的非零基底行元素 (phase-2)
			if (runMinSimplexMethod()) handleBound(); // 處理有界情況, 如果無界會提前終止並處理
		}
	}
	
	void print() { // debug
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
		
		printf("Type: "); // 印出 LP 的解的類型
		if (solutionType == Type::BOUNDED) printf("Bounded\n");
		else if (solutionType == Type::UNBOUNDED) printf("Unbounded\n");
		else if (solutionType == Type::INFEASIBLE) printf("Infeasible\n");
		
		printf(isMin ? "Minimum = " : "Maximum = "); // 印出極值
		printf("%.2f\n", extremum);
		
		printf("Solution: "); // 印出解向量
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
	VarBimap bimap;
	bool isMin; // min = 1 / max = 0
	Linearform objFunc; // 目標函數
	vector<Constraint> multiCon; // 多個約束

public:
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
	
	void test() { // test
		LP(isMin, objFunc, multiCon, {{ 0, FP64_INF }, { 0, FP64_INF }}).print();
	}
	
	void print() { // debug
		printf(isMin ? "min " : "max ");
		objFunc.print(bimap);
		for (Constraint& con: multiCon) con.print(bimap);
	}
};

class Tester { // 有一些範例用來測試正確性
private:
	void printTime(double t) { // print 一個時間 (輸入 ns)
		if (t < 1e3) {
			printf("%.f ns (1e-9 s)", t); return;
		};
		t /= 1e3;
		
		if (t < 1e3) {
			printf("%.f us (1e-6 s)", t); return;
		};
		t /= 1e3;
		
		if (t < 1e3) {
			printf("%.f ms (1e-3 s)", t); return;
		};
		t /= 1e3;
		
		printf("%.f s (1 s)", t);
	}
	
	struct TestLP { // LP 的測資結構
		bool isMin;
		Linearform objFunc;
		vector<Constraint> multiCon;
		vector<pair<double, double>> varRange;
	};
	
	vector<TestLP> lpTests = { // 所有 LP 的測資
		{ // test 0: Bounded
			false, Linearform().add(1, 0).add(1, 1), // max x0 + x1
			{
				Constraint().add(4, 0).add(3, 1).leq(17), // 4 x0 + 3 x1 <= 17
				Constraint().add(2, 0).add(-5, 1).geq(-9), // 2 x0 - 5 x1 >= -9
				Constraint().add(1, 0).add(10, 1).geq(25) // x0 + 10 x1 >= 25
			},
			{ { 0, FP64_INF }, { 0, FP64_INF } }
		},
		{ // test 1: Infeasible
			false, Linearform().add(1, 0).add(1, 1), // max x0 + x1
			{
				Constraint().add(4, 0).add(3, 1).leq(17), // 4 x0 + 3 x1 <= 17
				Constraint().add(2, 0).add(-5, 1).geq(-9), // 2 x0 - 5 x1 >= -9
				Constraint().add(1, 0).add(10, 1).geq(30) // x0 + 10 x1 >= 30
			},
			{ { 0, FP64_INF }, { 0, FP64_INF } }
		},
		{ // test 2: Unbounded (https://www.youtube.com/watch?v=113fYDGVzZk&t=30s)
			false, Linearform().add(1, 0), // max x0
			{
				Constraint().add(1, 0).add(-1, 1).leq(1), // x0 - x1 <= 1
				Constraint().add(2, 0).add(-1, 1).leq(4) // 2 x0 - x1 <= 4
			},
			{ { 0, FP64_INF }, { 0, FP64_INF } }
		},
		{ // test 3: Infeasible (edit by https://www.youtube.com/watch?v=58Jclpg0JCU&t=226s)
			false, Linearform().add(3, 0).add(1, 1), // max 3 x0 + x1
			{
				Constraint().add(4, 0).add(2, 1).leq(11), // 4 x0 + 2 x1 <= 11
			},
			{ { 2, FP64_INF }, { 2, FP64_INF } } // x0 >= 2, x1 >= 2
		},
		{ // test 4: Bounded (edit by https://www.youtube.com/watch?v=58Jclpg0JCU&t=226s)
			false, Linearform().add(3, 0).add(1, 1), // max 3 x0 + x1
			{
				Constraint().add(4, 0).add(2, 1).leq(11), // 4 x0 + 2 x1 <= 11
			},
			{ { 1, 2 }, { 0, FP64_INF } } // x0 <= [1, 2]
		}
	};

public:
	void testLp(uint32_t i) { // 測試 LP 的第 i 個側資, i = 0, 1, 2, ...
		if (i >= lpTests.size()) return; // 越界
		TestLP test = lpTests[i];
		
		auto start = chrono::high_resolution_clock::now(); // 測速
		LP lp(test.isMin, test.objFunc, test.multiCon, test.varRange); // 計算 LP 問題
		auto end = chrono::high_resolution_clock::now();
		
		printf("---------- test %d ----------\n", i);
		lp.print(); // print 一些資訊
		printf("Execution time: "); // LP 運算耗費的時間 (sec)
		printTime(chrono::duration<double, nano>(end - start).count());
		printf("\n\n");
	}
	
	void testAllLp() { // 測試所有 LP 的測資
		for (size_t i = 0; i < lpTests.size(); i++) testLp(i);
	}
};

int32_t main() {
	Tester tester;
	tester.testLp(4);
	
	/*
	IP("max", {{ 1, "x" }, { 1, "y" }})
		.addConstraint({{ 4, "x" }, { 3, "y" }}, "<=", 17)
		.addConstraint({{ 2, "x" }, { -5, "y" }}, ">=", -9)
		.addConstraint({{ 1, "x" }, { 10, "y" }}, ">=", 25).test();
	*/
	return 0;
}
