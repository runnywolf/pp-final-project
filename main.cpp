#include <stdio.h>
#include <cstdint>
#include <limits>
#include <vector>
#include <unordered_map>

using namespace std;

template <typename K, typename V>
bool mapHasKey(const unordered_map<K, V>& m, const K& key) {
	return m.find(key) != m.end(); // 如果 unordered_map 存在某個 key, 則回傳 value, 否則回傳 defaultValue
}

struct VarBimap { // 變數名稱與 x_j 的雙向映射
	unordered_map<string, uint32_t> strToIndex;
	unordered_map<uint32_t, string> indexToStr;
	
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
} bimap;

enum class Relation { LEQ, EQ, GEQ }; // <= (LEQ), = (EQ), >= (GEQ)

class Linearform { // 線性函數
private:
	unordered_map<uint32_t, double> form; // x_j index 映射到係數

public: // 有做 chaining
	Linearform& term(double coef, const string& varName) { // 添加一個未知數項, 變數名可以是任何字串
		if (coef != 0) form[bimap.getVarIndex(varName)] = coef; // 為字串變數分配一個 x_j index, 然後 map 到係數
		return *this;
	}
	
	void print() { // debug
		uint32_t c = 0;
		for (uint32_t i = 0; i < bimap.getVarCount(); i++) if (mapHasKey<uint32_t, double>(form, i)) {
			if (c++) printf(" + ");
			printf("%.2f[%s]", form[i], bimap.getVarName(i).c_str());
		}
	}
};

class Constraint: public Linearform { // 約束. 由一個線性函數, 關係, 和右側常數組成: c1 x1 + c2 x2 + ... ~ r
private:
	Relation relation; // >= or = or <=
	double rightConst; // 右側常數
	
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
	
	void print() { // debug
		Linearform::print();
		if (relation == Relation::LEQ) printf(" <= ");
		else if (relation == Relation::EQ) printf(" = ");
		else if (relation == Relation::GEQ) printf(" >= ");
		printf("%.2f", rightConst);
	}
};

class LinearProgramming { // LP
public:
	LinearProgramming(bool isMax, Linearform objectiveFunction) {
		
	}
};

int main() {
	numeric_limits<double>::infinity();
	
	Linearform().term(3, "x").term(0.1, "y").term(-2.33, "1z66A").print(); printf("\n");
	Constraint().term(3, "x").term(0, "y").term(-2.33, "z").geq(4.5).print(); printf("\n");
	
	return 0;
}
