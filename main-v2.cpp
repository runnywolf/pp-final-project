#include <stdio.h>
#include <cstdint>
#include <vector>
#include <unordered_map>

using namespace std;

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
	unordered_map<uint32_t, double> form; // x_j index 映射到係數
	
	Linearform& add(double coef, uint32_t varIndex) { // 添加一個未知數項 (chaining)
		if (mapHasKey<uint32_t, double>(form, varIndex)) form[varIndex] += coef; // x_j 存在就加上係數
		else form[varIndex] = coef; // x_j 不存在就建立未知數項
		
		return *this;
	}
	
	void negate() { // 對這個線性函數 *-1
		for (auto& term: form) term.second = -term.second; // 對所有係數變號
	}
	
	void print(VarBimap& bimap, bool addNewLine = true) { // debug
		uint32_t c = 0;
		for (uint32_t i = 0; i < bimap.getVarCount(); i++) if (mapHasKey<uint32_t, double>(form, i)) {
			if (c++) printf(" + ");
			printf("%.2f[%s]", form[i], bimap.getVarName(i).c_str());
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
	
	unordered_map<uint32_t, double> getLinearform() {
		return linearform.form;
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
	
	public:
		uint32_t rows; // row 數
		uint32_t cols; // col 數
		vector<int32_t> baseVarIndexs; // 基底變數的編號, -1 代表 slack var, -2 代表 artifical var
		
		double& operator()(uint32_t i, uint32_t j) { // 訪問扁平化的二維陣列
			return arr[cols * i + j];
		};
		
		void init(uint32_t rows, uint32_t cols) { // 初始化 tableau
			this->rows = rows; // 設定列數
			this->cols = cols; // 設定行數
			arr = vector<double>(rows * cols, 0); // 預設全為零
			baseVarIndexs = vector<int32_t>(cols, 0); // 因為第零列沒有基底編號, 使 index 對齊
		}
		
		void print() { // debug
			for (uint32_t i = 0; i < rows; i++) {
				for (uint32_t j = 0; j < cols; j++) printf("%.02f ", (*this)(i, j));
				if (i > 0) printf("= x_%d", baseVarIndexs[i]); // print 基底
				printf("\n");
			}
		}
	} tableau;
	
	void setTableauRow(uint32_t rowIndex, Constraint con, uint32_t& slackVarColIndex) { // 將一個約束加入到 tableau, 注意: slackVarColIndex 會被修改
		for (auto& [varIndex, coef]: con.getLinearform()) tableau(rowIndex, varIndex) = coef; // 填入一般變數
		
		if (con.haveSlackVar()) tableau(rowIndex, slackVarColIndex++) = con.getSlackVarCoef(); // 如果有 slack var, 需要添加係數 1 或 -1 到 tableau 內
		
		tableau(rowIndex, tableau.cols - 1) = con.getRightConst(); // 在列的最右元素, 設定右側常數 (基底的值)
		tableau.baseVarIndexs[rowIndex] = con.haveArtificialVar() ? -2 : -1; // -1 代表 slack var, -2 代表 artifical var
	}

public:
	LP(bool isMin, Linearform& objFunc, vector<Constraint>& multiCon, uint32_t varCount) {
		// TODO: 變數範圍應該傳入一個 tree node, 通過回溯祖先得到一串變數範圍 
		vector<Constraint> varRangeMultiCon = { Constraint().add(1, 0).geq(1) }; // 將變數範圍轉為多個約束
		
		// init tableau
		uint32_t slackVarCount = 0; // 計算 slack var 個數, 決定運算矩陣的 col 數 (因為要分配連續記憶體)
		for (Constraint& con: multiCon) if (con.haveSlackVar()) slackVarCount++;
		for (Constraint& con: varRangeMultiCon) if (con.haveSlackVar()) slackVarCount++;
		tableau.init(1 + multiCon.size() + varRangeMultiCon.size(), varCount + slackVarCount + 1); // init tableau
		
		// insert constraint into tableau
		uint32_t rowIndex = 1;
		uint32_t slackVarColIndex = varCount; // 因為一般變數的 col index 為 0 ~ varCount-1
		for (Constraint& con: multiCon) setTableauRow(rowIndex++, con, slackVarColIndex);
		for (Constraint& con: varRangeMultiCon) setTableauRow(rowIndex++, con, slackVarColIndex);
		
		// handling infeasible LP
		
		
		tableau.print();
	}
};

class IP { // Integer Programming
private:
	VarBimap bimap;
	bool isMin; // min = 1 / max = 0
	Linearform objFunc; // 目標函數
	vector<Constraint> multiCon; // 多個約束
	// TODO: var range tree

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
		LP(isMin, objFunc, multiCon, bimap.getVarCount());
	}
	
	void print() { // debug
		printf(isMin ? "min " : "max ");
		objFunc.print(bimap);
		for (Constraint& con: multiCon) con.print(bimap);
	}
};

int32_t main() {
	IP("max", {{ 1, "x" }, { 1, "y" }})
		.addConstraint({{ 4, "x" }, { 3, "y" }}, "<=", 17)
		.addConstraint({{ 2, "x" }, { -5, "y" }}, ">=", -9)
		.addConstraint({{ 1, "x" }, { 10, "y" }}, ">=", 25).test();
	
	return 0;
}
