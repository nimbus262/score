SCORER-Gap: Sequentially Correlated Rules forEvent Recommendation Considering Gap Size
===
This project contains the implementation and test setting of the *SCORER-Gap: Sequentially Correlated Rules for Event Recommendation Considering Gap Size*. On sparse datasets it can be assumed that rules occuring frequently with a low gap on average between the antecedent and consequent correlate together.
With SCORER-Gap we introduce an exponential weighting on the gapsize between the antecedent and the consequent of an occuring rule. An application on three real-world datasets shows that the mining procedure results in a different ruleset with slightly higher accuracy values.

Usage
===
* To see first results, run the script with the prepared parameters. These are modified to cope with the industrial dataset.
* The method "run_with_hyperparams()" can be generate different resulting rulesets based on a range of minsup, minconf and corrfac values at once.
* Plots are created to show the shift in the amount of rules and memory usage.
* Furthermore, in line 324f. the results can be compared regarding accuracy and coverage based on the strength value as an evaluation method.

References
===
* For further analysis of different real-world datasets, download the desired sequential database and adjust PATH and parameter variables.
* Besides the given datasets, please feel encouraged to also apply the implementation on real-world logs which you can find on the following site: [SPMF](https://www.philippe-fournier-viger.com/spmf/index.php?link=datasets.php).

Because of the adaption of an algorithm from the [SPMF framework](https://www.philippe-fournier-viger.com/spmf/index.php) this code is distributed under the GPL v3 license.