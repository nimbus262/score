#!/usr/bin/env python
# coding: utf-8


import os
import subprocess
import time
import random
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import pandas as pd

import srm_eval

PATH = "ca.pfv.spmf.algorithms.sequential_rules.rulegrowth"
# Data

TRAIN_KOSARAK = "../../../../database/kosarak_train.txt"
TRAIN_INDUST = "../../../../database/indust_train.txt"
TRAIN_COVID = "../../../../database/covid19_train.txt"

def get_current_dir():
    returned_output = subprocess.check_output("cd", shell=True,stderr=subprocess.STDOUT)
    print("Currently here: ", returned_output.decode("utf-8"))


results = []

# Test algorithm with different hyperparameter settings
def call(algo, sup, conf, data, out, consequent, corrfac):
    command = "java " + PATH + "." + algo + " {} {} {} {} {} {}".format(sup,conf,data,out,consequent,corrfac)
    try:
        returned_output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
    out_string = returned_output.decode("utf-8")
    print(out_string)
    mem = None
    nr_rules = None
    time = None
    for line in out_string.split("\n"):
        line = line.strip()
        if line.startswith("Sequential rules count:"):
            nr_rules = line.split(" ")[-1]
        if line.startswith("Total time:"):
            time = line.split(" ")[-2]
        if line.startswith("Max memory:"):
            mem = line.split(" ")[-1]
    results.append({"algorithm": algo, "support": sup, "confidence": conf, "corrfac": corrfac, 
                    "number of rules": nr_rules, "time": time, "memory": mem})
    with open("stats_" + out, "w") as text_file:
        text_file.write(returned_output.decode("utf-8"))
        
def run_with_hyperparams(algo, data, hyp_sup, hyp_conf, hyp_dec):
    print("Running " + algo)
    print(hyp_sup)
    print(hyp_conf)
    print(hyp_dec)
    counter = hyp_conf["amount"]*hyp_sup["amount"]*hyp_dec["amount"]
    print("Running " + str(counter) + " times with different parameters")
    for min_conf in np.linspace(hyp_conf["range"][0], hyp_conf["range"][1], hyp_conf["amount"]):
        for min_sup in np.linspace(hyp_sup["range"][0], hyp_sup["range"][1], hyp_sup["amount"]):
            for corrfac in np.linspace(hyp_dec["range"][0], hyp_dec["range"][1], hyp_dec["amount"]):
                start_time = time.time()
                print("### Current params: Minsup: {}, Minconf: {}, Corrfac: {}".format(min_sup, min_conf, corrfac))
                print(str(counter) + " computation(s) to go")
                out_label = data.split("/")[-1].split(".")[0] + "_" + algo + "_SUP{}_CONF{}_DECAY{}.txt".format(min_sup,min_conf,corrfac)
                out_file = "result/" + out_label
                #get_ipython().run_line_magic('cd', 'spmf/out/production/spmf')
                os.chdir('spmf/out/production/spmf')
                call(algo, min_sup, min_conf, data, out_label, 1, corrfac)
                #get_ipython().run_line_magic('cd', '../../../..')
                os.chdir('../../../..')
                #print(get_current_dir())
                if os.path.isfile(out_file): # if file exists, remove it
                    os.remove(out_file)
                if os.path.isfile("stats_" + out_file): # if file exists, remove it
                    os.remove("stats_" + out_file)
                os.rename("spmf/out/production/spmf/stats_" + out_label, "result/stats_" + out_label)
                os.rename("spmf/out/production/spmf/" + out_label, out_file)
                counter -= 1
                print("Computation took {} seconds".format(time.time()-start_time))
                
def filter_results(algo, conf, sup, x_val, y_val):
    x = []
    y = []
    for result in results:
        if result["algorithm"] == algo and result["confidence"] == conf and result["support"] == sup:
            x.append(float(result[x_val]))
            y.append(float(result[y_val]))
    return x, y


# In[ ]:

### The following creates all combinations of hyperparameter - use carefully
"""
run_with_hyperparams(algo="AlgoRULEGROWTH", data=TRAIN_COVID,  
                    hyp_sup={"range": [.2,.8], "amount": 4}, 
                    hyp_conf={"range": [.8,.8], "amount": 1},
                    hyp_dec={"range": [1,.2], "amount": 5}) 




run_with_hyperparams(algo="AlgoRULEGROWTH", data=TRAIN_KOSARAK,  
                    hyp_sup={"range": [.006,.012], "amount": 4}, 
                    hyp_conf={"range": [.8,.8], "amount": 1},
                    hyp_dec={"range": [1,.2], "amount": 5}) 


"""

run_with_hyperparams(algo="AlgoRULEGROWTH", data=TRAIN_INDUST,  
                    hyp_sup={"range": [.01,.01], "amount": 1}, 
                    hyp_conf={"range": [.7,.7], "amount": 1},
                    hyp_dec={"range": [.1,.1], "amount": 1}) 



def save_meta_results(results, dataset_title):
    with open("meta_results_" + dataset_title, "w") as text_file:
        for result in results:
            text_file.write(str(result)+"\n")

# Save results of current run        
#save_meta_results(results, "covid.txt")



# In[ ]:


# Plot number of rules and memory usage for different support values
"""
# Covid #rules
conf = .8

x, y = filter_results("AlgoRULEGROWTH", conf, .2, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .2", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .4, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .4", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, 0.6000000000000001, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .6", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .8, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .8", marker="p")
plt.title("Amount of rules with RuleGrowth and Correlation Gap Constraint (Conf=.8)")
plt.xlabel('Correlation Factor')
plt.ylabel('Number of Rules')
plt.legend()
plt.show()
plt.close()


# In[ ]:


# Covid memory usage
conf = .8
axes = plt.gca()
axes.set_ylim([20,22])
x, y = filter_results("AlgoRULEGROWTH", conf, .2, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .2", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .4, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .4", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, 0.6000000000000001, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .6", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .8, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .8", marker="p")
plt.title("Memory results with RuleGrowth and Correlation Gap Constraint (Conf=.8)", pad=20)
plt.xlabel('Correlation Factor')
plt.ylabel('Memory (MB)')
plt.legend()
plt.show()
plt.close()


# In[ ]:


# Kosarak #rules
conf = .8

x, y = filter_results("AlgoRULEGROWTH", conf, .006, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .006", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .008, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .008", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, 0.01, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .01", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .012, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .012", marker="p")
plt.title("Amount of rules with RuleGrowth and Correlation Gap Constraint (Conf=.8)")
plt.xlabel('Correlation Factor')
plt.ylabel('Number of Rules')
plt.legend()
plt.show()
plt.close()


# In[ ]:


# Kosarak memory usage
conf = .8

x, y = filter_results("AlgoRULEGROWTH", conf, .006, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .006", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .008, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .008", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, .01, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .01", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .012, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .012", marker="p")
plt.title("Memory results with RuleGrowth and Correlation Gap Constraint (Conf=.8)", pad=20)
plt.xlabel('Correlation Factor')
plt.ylabel('Memory (MB)')
plt.legend()
plt.show()
plt.close()


# In[ ]:


# Industrial dataset #rules
conf = .8

x, y = filter_results("AlgoRULEGROWTH", conf, .002, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .002", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .004, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .004", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, 0.006, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .006", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .008, "corrfac", "number of rules")
plt.plot(x, y, label = "Minimum Support = .008", marker="p")
plt.title("Amount of rules with RuleGrowth and Correlation Gap Constraint (Conf=.8)")
plt.xlabel('Correlation Factor')
plt.ylabel('Number of Rules')
plt.legend()
plt.show()
plt.close()


# In[ ]:


# Industrial dataset memory usage
conf = .8

x, y = filter_results("AlgoRULEGROWTH", conf, .002, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .002", marker="o")
x, y = filter_results("AlgoRULEGROWTH", conf, .004, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .004", marker="^")
x, y = filter_results("AlgoRULEGROWTH", conf, .006, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .006", marker="s")
x, y = filter_results("AlgoRULEGROWTH", conf, .008, "corrfac", "memory")
plt.plot(x, y, label = "Minimum Support = .008", marker="p")
plt.title("Memory results with RuleGrowth and Correlation Gap Constraint (Conf=.8)", pad=20)
plt.xlabel('Correlation Factor')
plt.ylabel('Memory (MB)')
plt.legend()
plt.show()
plt.close()
"""

# # Analyze Output


def evaluate(outfile, test_path, eval_method="confidence", gui=False):
    evaluation = srm_eval.SRMEval("result/" + outfile, gui)
    df = evaluation.to_dataframe(gui)
    coverage = None
    accuracy = None
    # Test performance according to evaluation metric
    if eval_method == "strength":
        coverage, accuracy = srm_eval.test_performance(df, test_path, "strength", gui)
    else:
        coverage, accuracy = srm_eval.test_performance(df, test_path, "confidence", gui)
    return df, accuracy, coverage
    

def filter_rules_for_recommendation(df, corrfac, eval_method="strength", gui=False):
    conf = list(df["conf"])
    if not gui:
        decsup = list(df["decsup"])
    sup = list(df["sup"])
    rules = list(df["rule"])
    #print(df)
    strength = []
    for index, row in df.iterrows():
    #    print(row['sup'], row['conf'])
        if not corrfac:
            strength.append(row['sup'] * row['conf'])
        else:
            strength.append(row['decsup'] * row['conf'])
        
    #print(strength)
    df['strength'] = np.array(strength)
    result = None
    if eval_method == "strength":
        result = df.sort_values("strength", ascending=False)
    elif eval_method == "confidence":
        result = df.sort_values("conf", ascending=False)
    elif eval_method == "support":
        result = df.sort_values("sup", ascending=False)
    #print(df)
    return result

def parse_rule(rule_string):
    antecedent = []
    consequent = []
    rule_arr = rule_string.strip().split("==>")
    antecedent = rule_arr[0]
    consequent = rule_arr[-1]
    return [antecedent, consequent]

def equals(r1, r2):
    return r1[0] == r2[0] and r1[1] == r2[1]


# Comparison of SCORER-Gap and Original Approach

PATH_RESULT1 = "indust_train_AlgoRULEGROWTH_SUP0.01_CONF0.7_DECAY0.1.txt"
PATH_RESULT2 = "indust_train_AlgoRULEGROWTH_SUP0.01_CONF0.7_DECAY0.1.txt"

print("Compare results")
df1, accuracy1, coverage1 = evaluate(PATH_RESULT1, 
                                  "database/indust_test.txt",
                                  eval_method="strength", gui=False)
print(df1)
sorted_df1 = filter_rules_for_recommendation(df1, corrfac=False)


df3, accuracy3, coverage3 = evaluate(PATH_RESULT2, 
                                  "database/indust_test.txt",
                                  eval_method="strength", gui=False)
print(df3)
sorted_df3 = filter_rules_for_recommendation(df3, corrfac=False)

comp_n = 20
print("Comparing first {} rows:\n{}\n{}".format(comp_n, sorted_df1.head(comp_n)[["rule", "strength"]], sorted_df3.head(comp_n)[["rule", "strength"]]))


# In[25]:

ruleset_corrfac = set(sorted_df1.head(58)["rule"].values)
ruleset_orig = set(sorted_df3.head(58)["rule"].values)
different_rules = ruleset_corrfac.symmetric_difference(ruleset_orig)
unique_corrfac = ruleset_corrfac.difference(ruleset_orig)
unique_orig = ruleset_orig.difference(ruleset_corrfac)
common_rules = ruleset_corrfac.intersection(ruleset_orig)
print(len(common_rules))
#print(common_rules)
print(len(unique_corrfac))
print(len(unique_orig))
print(len(different_rules))
#print(unique_corrfac)


# In[ ]:

# Write common and different rules to file
with open("common_rules.txt", "w+") as cr:
    for rule in common_rules:
        row = sorted_df1.loc[sorted_df1['rule'] == rule]
        #print(row["sup"].values[0])
        cr.write(rule + ", SUP: " + str(row["sup"].values[0]) + ", CONF: " + str(row["conf"].values[0]) + ", CORRSUP: " + str(row["decsup"].values[0]) + ", STR: " + str(row["strength"].values[0]) + "\n")
with open("unique_corrfac.txt", "w+") as cr:
    for rule in unique_corrfac:
        row = sorted_df1.loc[sorted_df1['rule'] == rule]
        cr.write(rule + ", SUP: " + str(row["sup"].values[0]) + ", CONF: " + str(row["conf"].values[0]) + ", CORRSUP: " + str(row["decsup"].values[0]) + ", STR: " + str(row["strength"].values[0]) + "\n")
with open("unique_orig.txt", "w+") as cr:
    for rule in unique_orig:
        row = sorted_df3.loc[sorted_df3['rule'] == rule]
        cr.write(rule + ", SUP: " + str(row["sup"].values[0]) + ", CONF: " + str(row["conf"].values[0]) + ", CORRSUP: " + str(row["decsup"].values[0]) + ", STR: " + str(row["strength"].values[0]) + "\n")


# In[ ]:


print("Stats of rules which both approaches have in common in the result set\n")
common_rules_df = sorted_df1[(sorted_df1['rule'].isin(common_rules))]
#print(common_rules_df)
#print(common_rules_df[['sup', 'conf', 'decsup', 'strength']].mean())
print(common_rules_df.describe())
print("##################\n")
print("Stats of rules which are only in the result set of scorer-gap approach\n")
unique_corrfac_df = sorted_df1[(sorted_df1['rule'].isin(unique_corrfac))]
#print(unique_corrfac_df[['sup', 'conf', 'decsup', 'strength']].mean())
print(unique_corrfac_df.describe())
print("##################\n")
print("Stats of rules which are only in the result set of rulegrowth approach\n")
unique_orig_df = sorted_df3[(sorted_df3['rule'].isin(unique_orig))]
#print(unique_orig_df[['sup', 'conf', 'decsup', 'strength']].mean())
print(unique_orig_df.describe())


# In[ ]:


mean_antecedent_len_orig = unique_orig_df.antecedent.str.len().sum()/len(unique_corrfac_df.index)
mean_antecedent_len_corrfac = unique_corrfac_df.antecedent.str.len().sum()/len(unique_corrfac_df.index)
print("Mean length of antecedent in the rulegrowth result set")
print(mean_antecedent_len_orig)
print("Mean length of antecedent in the scorer-gap result set")
print(mean_antecedent_len_corrfac)


# Analyze Database

# create list of sequences (as list) from input text file
def preprocess_db(path, sample=False, sample_size=.3):
    new_db = []
    with open(path) as file:
        db = file.readlines()
        print("Example sequence:\n" + db[0])
        if sample:
            db = get_sample(db, sample_size)
        for sequence in db:
            new_sequence= []
            sequence = sequence.split()
            for item in sequence:
                item = int(item)
                #print(item)
                if item != -1 and item != -2:
                    new_sequence.append(item)
            new_db.append(new_sequence)
    return new_db

# sample_size is a floating point number representing the mount of the original database (in percent)
def get_sample(db, sample_size):
    sample = []
    for sequence in db:
        if random.random() < sample_size:
            sample.append(sequence)
    #print(len(sample))
    return sample
            
# with std
def get_nr_of_seq(db):
    return len(db)

def get_nr_of_dist_items(db):
    flat_list = [item for sublist in db for item in sublist]
    db_as_set = set(flat_list)
    #print(db_as_set)
    return len(db_as_set)

# with std
def get_mean_seq_length(db):
    lengths = []
    for sequence in db:
        lengths.append(len(sequence))
    return sum(lengths)/len(db), np.std(lengths)

# with std
def get_nr_of_dist_items_per_seq(db):
    dist_items = []
    for sequence in db:
        seq_set = set(sequence)
        dist_items.append(len(seq_set))
    return sum(dist_items)/len(db), np.std(dist_items)


# In[ ]:


PATH_KOSARAK = "database/db_kosarak.txt"
PATH_INDUST = "database/db_indust.txt"
PATH_COVID19 = "database/db_MT745584.txt"


# In[ ]:


print("Analyzing Covid19 data")
prep_db = preprocess_db(PATH_COVID19)
print("Number of sequences: ", get_nr_of_seq(prep_db))
print("Number of distinct items: ", get_nr_of_dist_items(prep_db))
print("Average sequence length: ", get_mean_seq_length(prep_db))
print("Number of distinct items per sequence: ", get_nr_of_dist_items_per_seq(prep_db))


# In[ ]:


print("Analyzing industrial data")
prep_db = preprocess_db(PATH_INDUST)
print("Number of sequences: ", get_nr_of_seq(prep_db))
print("Number of distinct items: ", get_nr_of_dist_items(prep_db))
print("Average sequence length: ", get_mean_seq_length(prep_db))
print("Number of distinct items per sequence: ", get_nr_of_dist_items_per_seq(prep_db))


# In[ ]:


print("Analyzing Kosarak data")
prep_db = preprocess_db(PATH_KOSARAK, sample=True, sample_size=.09)
print("Number of sequences: ", get_nr_of_seq(prep_db))
print("Number of distinct items: ", get_nr_of_dist_items(prep_db))
print("Average sequence length: ", get_mean_seq_length(prep_db))
print("Number of distinct items per sequence: ", get_nr_of_dist_items_per_seq(prep_db))
