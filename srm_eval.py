import pandas as pd
import time
import os

class SRMEval:
    def __init__(self, result_file, gui):
        self.rule_list = []
        self.antecedent_list = []
        self.consequent_list = []
        self.sup_list = []
        self.conf_list = []
        self.decsup_list = []
        self.meangap_list = []
        
        with open(result_file) as file:
            number_of_rules = 0
            for row in file:
                number_of_rules += 1
                if gui:
                    sharp_1 = 0
                    sharp_2 = 0
                    for i, char in enumerate(row):
                        if char =="#" and row[i+1]=="S": # support
                            sharp_1 = i
                        if char =="#" and row[i+1]=="C": # confidence
                            sharp_2 = i
                    rule = row[:sharp_1]
                    self.rule_list.append(rule)
                    splitted_r = rule.split(" ")
                    antecedent = splitted_r[0]
                    consequent = splitted_r[2]
                    self.antecedent_list.append([int(elem) for elem in antecedent.split(",")])
                    self.consequent_list.append(int(consequent))
                    sup = int(row[sharp_1+6:sharp_2])
                    self.sup_list.append(sup)
                    conf = float(row[sharp_2+7:])
                    self.conf_list.append(conf)
                else:
                    sharp_1 = 0
                    sharp_2 = 0
                    sharp_3 = 0
                    sharp_4 = 0
                    for i, char in enumerate(row):
                        if char =="#" and row[i+1]=="S": # support
                            sharp_1 = i
                        if char =="#" and row[i+1]=="D": # decayed support
                            sharp_2 = i
                        if char =="#" and row[i+1]=="C": # confidence
                            sharp_3 = i
                        if char =="#" and row[i+1]=="M": # confidence
                            sharp_4 = i
                    rule = row[:sharp_1]
                    self.rule_list.append(rule)
                    splitted_r = rule.split(" ")
                    antecedent = splitted_r[0]
                    consequent = splitted_r[2]
                    self.antecedent_list.append([int(elem) for elem in antecedent.split(",")])
                    self.consequent_list.append(int(consequent))
                    sup = int(row[sharp_1+6:sharp_2])
                    self.sup_list.append(sup)
                    decsup = float(row[sharp_2+9:sharp_3])
                    self.decsup_list.append(decsup)
                    conf = float(row[sharp_3+7:sharp_4])
                    self.conf_list.append(conf)
                    meangap = float(row[sharp_4+10:])
                    self.meangap_list.append(meangap)
            print("Got {} rules".format(number_of_rules))

    def to_dataframe(self, gui):
        # combine everything into a dataframe
        if gui:
            df = pd.DataFrame(list(zip(self.rule_list, self.antecedent_list, self.consequent_list, self.sup_list, self.conf_list)), 
                    columns =['rule', 'antecedent', 'consequent', 'sup', 'conf'])
        else:
            df = pd.DataFrame(list(zip(self.rule_list, self.antecedent_list, self.consequent_list, self.sup_list, self.conf_list, self.decsup_list, self.meangap_list)), 
                    columns =['rule', 'antecedent', 'consequent', 'sup', 'conf', 'decsup', 'meangap'])
        return df
    


def test_performance(df, test_data, eval_method="strength", gui=False):
    start_time = time.time()
    conf = list(df["conf"])
    if not gui:
        decsup = list(df["decsup"])
    sup = list(df["sup"])
    rules = list(df["rule"])
    consequent_list = list(df["consequent"])
    antecedent_list = list(df["antecedent"])
    # Read-in Test data
    test_list = read_test_data(test_data) 

    match_count = 0
    correct_pred = 0
    for test in test_list:
        matches = []
        for i in range(len(antecedent_list)):
            indices_antecedent = contains(antecedent_list[i], test)
            indices_consequent = contains(consequent_list[i], test)
            if -1 not in indices_antecedent:
                matches.append((consequent_list[i], sup[i], conf[i], sup[i] * conf[i], indices_antecedent, indices_consequent))
        if matches:
            # Different evaluation methods
            match_count +=1     
            if eval_method == "strength":
                matches.sort(key=lambda tup: tup[3])  # sorts in place
            elif eval_method == "confidence":
                matches.sort(key=lambda tup: tup[2])  # sorts in place
            elif eval_method == "support":
                matches.sort(key=lambda tup: tup[1])  # sorts in place
            else:
                raise Exception("Error: Could not understand evaluation method")
            best = matches[-1] # has been sorted ascending

            if best[0] in test and min(best[5]) > max(best[4]): # better use dictionary
                correct_pred += 1

    coverage = match_count/len(test_list)
    accuracy = correct_pred/match_count
    print("Calculation time: {}".format(time.time() - start_time))
    print("Metric: {}, Coverage: {}, Accuracy: {}\n".format(eval_method, coverage, accuracy))
    return (coverage, accuracy)


def contains(small, big):
    res = [] 
    new_small = []
    if not isinstance(small, list):
        new_small.append(small)
    else:
        new_small = small
    for el in new_small:
        indices = [i for i, x in enumerate(big) if x == el]
        if len(indices) == 0:
            res.append(-1)
        else:
            res.extend(indices)
    return res


def read_test_data(test_data):
    test_list = []    
    with open(test_data, 'r') as txtfile:
        for i, sequence in enumerate(txtfile):
            seq_str = ""
            for num in sequence.split(" "):
                num = num.strip()
                if str(num) == "-1" or str(num) == "-2\n":  
                    continue
                else:
                    if not num == "":
                        seq_str = seq_str + num + ","
            seq_str = seq_str[:-1]
            test_list.append([int(elem) for elem in seq_str.split(",")]) 
    return test_list


