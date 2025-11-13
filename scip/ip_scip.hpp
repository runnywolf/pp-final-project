#pragma once
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <iostream>

class IP {
private:
    SCIP* scip_;
    std::map<std::string, SCIP_VAR*> vars_;
    bool is_max_;
    
public:
    // 建構子: sense 可以是 "max" 或 "min"
    IP(const std::string& sense, 
       const std::vector<std::pair<double, std::string>>& obj) {
        
        // 初始化 SCIP
        SCIP_CALL_ABORT(SCIPcreate(&scip_));
        SCIP_CALL_ABORT(SCIPincludeDefaultPlugins(scip_));
        SCIP_CALL_ABORT(SCIPcreateProbBasic(scip_, "supply_chain"));
        
        // 設定目標方向
        is_max_ = (sense == "max" || sense == "maximize");
        SCIP_CALL_ABORT(SCIPsetObjsense(scip_, 
            is_max_ ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE));
        
        // 建立所有變數 (預設為非負整數)
        for (const auto& term : obj) {
            const std::string& var_name = term.second;
            if (vars_.find(var_name) == vars_.end()) {
                SCIP_VAR* var;
                SCIP_CALL_ABORT(SCIPcreateVarBasic(
                    scip_,
                    &var,
                    var_name.c_str(),
                    0.0,                    // lower bound
                    SCIPinfinity(scip_),    // upper bound (無限大)
                    term.first,             // objective coefficient
                    SCIP_VARTYPE_INTEGER    // 整數變數
                ));
                SCIP_CALL_ABORT(SCIPaddVar(scip_, var));
                vars_[var_name] = var;
            } else {
                // 變數已存在,累加目標係數
                SCIP_CALL_ABORT(SCIPaddVarObj(scip_, vars_[var_name], term.first));
            }
        }
    }
    
    // 解構子
    ~IP() {
        // 釋放所有變數
        for (auto& pair : vars_) {
            SCIP_CALL_ABORT(SCIPreleaseVar(scip_, &pair.second));
        }
        SCIP_CALL_ABORT(SCIPfree(&scip_));
    }
    
    // 新增限制式
    void addConstraint(const std::vector<std::pair<double, std::string>>& lhs,
                       const std::string& sense,
                       double rhs) {
        
        // 確保所有變數都存在
        for (const auto& term : lhs) {
            const std::string& var_name = term.second;
            if (vars_.find(var_name) == vars_.end()) {
                // 新增缺少的變數 (目標係數為 0)
                SCIP_VAR* var;
                SCIP_CALL_ABORT(SCIPcreateVarBasic(
                    scip_,
                    &var,
                    var_name.c_str(),
                    0.0,
                    SCIPinfinity(scip_),
                    0.0,  // 目標係數 = 0
                    SCIP_VARTYPE_INTEGER
                ));
                SCIP_CALL_ABORT(SCIPaddVar(scip_, var));
                vars_[var_name] = var;
            }
        }
        
        // 建立限制式
        SCIP_CONS* cons;
        
        if (sense == "=") {
            SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(
                scip_, &cons, "eq_cons",
                0, nullptr, nullptr,
                rhs, rhs  // lhs = rhs
            ));
        } else if (sense == "<=") {
            SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(
                scip_, &cons, "le_cons",
                0, nullptr, nullptr,
                -SCIPinfinity(scip_), rhs
            ));
        } else if (sense == ">=") {
            SCIP_CALL_ABORT(SCIPcreateConsBasicLinear(
                scip_, &cons, "ge_cons",
                0, nullptr, nullptr,
                rhs, SCIPinfinity(scip_)
            ));
        } else {
            throw std::runtime_error("Unknown sense: " + sense);
        }
        
        // 加入變數到限制式
        for (const auto& term : lhs) {
            SCIP_CALL_ABORT(SCIPaddCoefLinear(scip_, cons, 
                vars_[term.second], term.first));
        }
        
        SCIP_CALL_ABORT(SCIPaddCons(scip_, cons));
        SCIP_CALL_ABORT(SCIPreleaseCons(scip_, &cons));
    }
    
    // 求解
    bool solve() {
        SCIP_CALL_ABORT(SCIPsolve(scip_));
        
        SCIP_SOL* sol = SCIPgetBestSol(scip_);
        return (sol != nullptr);
    }
    
    // 取得目標值
    double getObjValue() const {
        return SCIPgetPrimalbound(scip_);
    }
    
    // 取得變數值
    double getVarValue(const std::string& var_name) const {
        auto it = vars_.find(var_name);
        if (it == vars_.end()) {
            return 0.0;
        }
        
        SCIP_SOL* sol = SCIPgetBestSol(scip_);
        if (sol == nullptr) {
            return 0.0;
        }
        
        return SCIPgetSolVal(scip_, sol, it->second);
    }
    
    // 取得求解狀態
    std::string getStatus() const {
        SCIP_STATUS status = SCIPgetStatus(scip_);
        switch(status) {
            case SCIP_STATUS_OPTIMAL: return "optimal";
            case SCIP_STATUS_INFEASIBLE: return "infeasible";
            case SCIP_STATUS_UNBOUNDED: return "unbounded";
            case SCIP_STATUS_TIMELIMIT: return "timelimit";
            default: return "unknown";
        }
    }
    
    // 設定時間限制 (秒)
    void setTimeLimit(double seconds) {
        SCIP_CALL_ABORT(SCIPsetRealParam(scip_, "limits/time", seconds));
    }
    
    // 設定輸出等級 (0 = 安靜, 4 = 詳細)
    void setVerbosity(int level) {
        SCIP_CALL_ABORT(SCIPsetIntParam(scip_, "display/verblevel", level));
    }
    
    // 印出所有變數值
    void printSolution() const {
        SCIP_SOL* sol = SCIPgetBestSol(scip_);
        if (sol == nullptr) {
            std::cout << "No solution found" << std::endl;
            return;
        }
        
        std::cout << "Objective value: " << getObjValue() << std::endl;
        std::cout << "Variables:" << std::endl;
        
        for (const auto& pair : vars_) {
            double val = SCIPgetSolVal(scip_, sol, pair.second);
            if (val > 1e-6) {  // 只印非零值
                std::cout << "  " << pair.first << " = " << val << std::endl;
            }
        }
    }
};