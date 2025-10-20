#include <stdio.h>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <unordered_map>

using namespace std;

enum class Relation { LEQ, EQ, GEQ }; // <= (LEQ), = (EQ), >= (GEQ)

template <typename K, typename V>
bool mapHasKey(const unordered_map<K, V>& m, const K& key) {
	return m.find(key) != m.end(); // 如果 unordered_map 存在某個 key, 則回傳 value, 否則回傳 defaultValue
}

struct VarBimap { // 變數名稱與 x_j 的雙向映射 (全域, 唯一)
	unordered_map<string, uint32_t> strToIndex;
	unordered_map<uint32_t, string> indexToStr;
	uint32_t slackVarCount = 0; // slack var 個數
	uint32_t artificialVarCount = 0; // artificial var 個數
	
	uint32_t getVarCount() { // 獲取已註冊變數的數量
		return strToIndex.size();
	}
	
	uint32_t getVarIndex(const string& varName) { // 註冊一個字串變數, 自動分配編號, 回傳這個字串變數對應的編號 j (x_j)
		if (mapHasKey<string, uint32_t>(strToIndex, varName)) return strToIndex[varName]; // 如果字串變數已註冊過, 回傳對應的編號 j
		
		const uint32_t newVarIndex = getVarCount(); // 如果字串變數沒有註冊過, 分配編號 0, 1, 2, ...
		strToIndex[varName] = newVarIndex; // 註冊雙向映射
		indexToStr[newVarIndex] = varName;
		return newVarIndex; // 回傳分配的新編號
	}
	
	string getVarName(uint32_t varIndex) { // 編號 j (x_j) 轉字串變數
		if (mapHasKey<uint32_t, string>(indexToStr, varIndex)) return indexToStr[varIndex];
		return "[unknown-var]";
	}
	
	string createSlackVar() { // 註冊一個 slack var, 並回傳 name & index
		const string varName = "@slk-" + to_string(++slackVarCount); // 產生一個獨立的 slack var name
		const uint32_t varIndex = getVarIndex(varName); // 註冊, 並獲取 var index
		return varName;
	}
	
	string createArtificialVar() { // 註冊一個 artificial var, 並回傳 name & index
		const string varName = "@art-" + to_string(++artificialVarCount); // 產生一個獨立的 slack var name
		const uint32_t varIndex = getVarIndex(varName); // 註冊, 並獲取 var index
		return varName;
	}
} bimap;

class Linearform { // 線性函數
protected:
	unordered_map<uint32_t, double> form; // x_j index 映射到係數

public: // 有做 chaining
	Linearform& term(double coef, const string& varName) { // 添加一個未知數項, 變數名可以是任何字串
		if (coef != 0) form[bimap.getVarIndex(varName)] = coef; // 為字串變數分配一個 x_j index, 然後 map 到係數
		return *this;
	}
	
	void negate() { // 對這個線性函數 *-1
		for (auto& term: form) term.second = -term.second; // 對所有係數變號
	}
	
	void print(bool addNewLine = true) { // debug
		uint32_t c = 0;
		for (uint32_t i = 0; i < bimap.getVarCount(); i++) if (mapHasKey<uint32_t, double>(form, i)) {
			if (c++) printf(" + ");
			printf("%.2f[%s]", form[i], bimap.getVarName(i).c_str());
		}
		if (addNewLine) printf("\n");
	}
};

class Constraint: public Linearform { // 約束. 由一個線性函數, 關係, 和右側常數組成: c1 x1 + c2 x2 + ... ~ r
private:
	double rightConst; // 右側常數
	Relation relation; // >= or = or <=
	
	void set(Relation relation, double rightConst) { // leq, eq, geq 共用的 set 邏輯
		this->relation = relation;
		this->rightConst = rightConst;
	}

public: // 有做 chaining
	Constraint& term(double coef, const string& varName) { // 添加一個未知數項, 變數名可以是任何字串
		Linearform::term(coef, varName); // 呼叫父類的版本
		return *this;
	}
	
	Constraint& leq(double rightConst) { // 添加最右側的 "<= r" 部分
		set(Relation::LEQ, rightConst);
		return *this;
	}
	
	Constraint& eq(double rightConst) { // 添加最右側的 "= r" 部分
		set(Relation::EQ, rightConst);
		return *this;
	}
	
	Constraint& geq(double rightConst) { // 添加最右側的 ">= r" 部分
		set(Relation::GEQ, rightConst);
		return *this;
	}
	
	void stdOfNegativeRightConst() { // 對負的右側常數進行標準化
		if (rightConst >= 0) return; // 若右側常數 >= 0, 跳過這一步
		
		rightConst = -rightConst; // 對右側常數變號
		negate(); // 對所有係數變號
		
		if (relation == Relation::LEQ) relation = Relation::GEQ; // 兩側變號需要將 <= 和 >= 轉向
		else if (relation == Relation::GEQ) relation = Relation::LEQ;
	}
	
	int32_t stdOfInequality() { // 對不等式進行標準化, 若 slack var 係數 = 1 則回傳 var index, 否則回傳 -1 (代表要添加 arti-var)
		if (relation == Relation::EQ) return -1; // "=" 不用添加 slack var
		
		const string slackVarName = bimap.createSlackVar(); // 產生一個獨立的 slack var name
		if (relation == Relation::LEQ) term(1, slackVarName); // "... <= r" -> "... + x = r"
		if (relation == Relation::GEQ) term(-1, slackVarName); // "... >= r" -> "... - x = r"
		
		const int32_t slackVarIndex = (relation == Relation::LEQ) ? bimap.getVarIndex(slackVarName) : -1; // 只有 <= 的 slack var 係數 = 1
		relation = Relation::EQ; // 這一步只是形式上的, 不影響後續演算法的正確性
		return slackVarIndex;
	}
	
	void print(bool addNewLine = true) { // debug
		Linearform::print(false);
		if (relation == Relation::LEQ) printf(" <= ");
		else if (relation == Relation::EQ) printf(" = ");
		else if (relation == Relation::GEQ) printf(" >= ");
		printf("%.2f", rightConst);
		if (addNewLine) printf("\n");
	}
};

class LinearProgramming { // LP
private:
	bool isMax; // max -> true, min -> false
	bool isFlip = false; // 如果 min 有被翻轉成 max 過
	Linearform objectiveFunction; // 目標函數
	vector<Constraint> multiCon; // 多個約束
	vector<int32_t> baseVarIndexs; // 基底, 多個 x_j = ... 的編號 j, 在執行 .std 之前是空的
	
	void stdOfMinMax() { // 透過變號操作將 min LP 轉為 max LP, isFlip 會保留 min/max 資訊
		if (!isMax) {
			objectiveFunction.negate();
			isMax = true;
			isFlip = true;
		}
	}
	
	void stdOfInfeasible() { // slack var 不足以形成基底, 需要添加額外的 artificial var
		Linearform minObjFunc; // 因為要找出一個最小的起始向量, 這邊會從求 max 改為求 min
		for (size_t i = 0; i < baseVarIndexs.size(); i++) if (baseVarIndexs[i] == -1) { // 如果編號為 -1 代表要添加 artificial var
			const string artificialVarName = bimap.createArtificialVar(); // 產生一個獨立的 artificial var name
			multiCon[i].term(1, artificialVarName); // 添加一個 artificial var 作為暫時性基底
			minObjFunc.term(1, artificialVarName); // 在 minimum function 也添加 artificial var
			baseVarIndexs[i] = bimap.getVarIndex(artificialVarName); // 紀錄 artificial var 的基底
		}
		
		if (bimap.artificialVarCount == 0) return; // 若 slack var 足夠, 跳過這一步
		
		// 轉 array
	}

public:
	LinearProgramming(bool isMax, Linearform& objectiveFunction) { // 宣告 LP 問題的 max/min 和 目標函數
		this->isMax = isMax;
		this->objectiveFunction = objectiveFunction;
	}
	
	LinearProgramming& addCon(Constraint& con) { // 添加約束 (chaining)
		multiCon.push_back(con);
		return *this;
	}
	
	void std() { // 開始執行標準化, 不等式會被修改
		for (Constraint& con: multiCon) {
			con.stdOfNegativeRightConst(); // 對負的右側常數進行標準化
			baseVarIndexs.push_back(con.stdOfInequality()); // 不存在 slack var 的 con index (改為複寫 vector, 因為要記錄右下向量)
		}
		stdOfMinMax(); // 透過變號操作將 min LP 轉為 max LP, isFlip 會保留 min/max 資訊
		stdOfInfeasible(); // slack var 不足以形成基底, 需要添加額外的 artificial var
	}
	
	void print() { // debug
		printf(isMax ? "max " : "min ");
		objectiveFunction.print();
		for (Constraint& con: multiCon) con.print();
	}
};

int32_t main() {
	LinearProgramming lp = LinearProgramming(false, Linearform().term(1, "x").term(1, "y"))
		.addCon(Constraint().term(4, "x").term(3, "y").leq(17))
		.addCon(Constraint().term(2, "x").term(-5, "y").geq(-9))
		.addCon(Constraint().term(1, "x").term(10, "y").geq(25));
	
	lp.print();
	lp.std();
	lp.print();
	
	return 0;
}
