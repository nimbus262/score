package ca.pfv.spmf.algorithms.sequential_rules.rulegrowth;
/* This file is copyright (c) 2008-2013 Philippe Fournier-Viger
 *
 * This file is part of the SPMF DATA MINING SOFTWARE
 * (http://www.philippe-fournier-viger.com/spmf).
 *
 * SPMF is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with
 * SPMF. If not, see <http://www.gnu.org/licenses/>.
 */

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

import ca.pfv.spmf.input.sequence_database_list_integers.Sequence;
import ca.pfv.spmf.input.sequence_database_list_integers.SequenceDatabase;
import ca.pfv.spmf.tools.MemoryLogger;

/**
 * This is the original implementation of the RULEGROWTH algorithm for mining sequential rules
 * common to several sequences where antecedent and consequent are unordered itemsets. The RuleGrowth
 * algorithm is described in this paper:
 * <br/><br/>
 *  Fournier-Viger, P., Nkambou, R. & Tseng, V. S. (2011).
 *  RuleGrowth: Mining Sequential Rules Common to Several Sequences by Pattern-Growth.
 *  Proceedings of the 26th Symposium on Applied Computing (ACM SAC 2011). ACM Press, pp. 954-959.
 * <br/><br/>
 * The main method of this algorithm is "runAlgorithm". It output the result to a file.
 *
 * @see Occurence
 * @see Sequence
 * @see SequenceDatabase
 * @author Philippe Fournier-Viger
 */
public class AlgoRULEGROWTH {
	//*** for statistics ***/

	/** start time of latest execution */
	long timeStart = 0;

	/**  end time of latest execution */
	long timeEnd = 0;

	/** number of rules generated */
	int ruleCount;

	//*** parameters ***/
	/** minimum confidence */
	double minConfidence;

	/** minimum absolute support */
	double minsupAbsolute;

	/** minimum relative support */
	double minsupRelative;

	/** this is the sequence database */
	SequenceDatabase database;

	/*** internal variables
	 // This map contains for each item (key) a map of occurences (value).
	 // The map of occurences associates to sequence ID (key), an occurence of the item (value). */
	Map<Integer,  Map<Integer, Occurence>> mapItemCount;  // item, <tid, occurence>

	/** object to write the output file */
	BufferedWriter writer = null;

	/** Used by the debug mode to keep all rules found */
	static List<Rule> allRulesFoundForDEBUG = new ArrayList<Rule>();

	/** If true the debug mode will be used */
	boolean DEBUG = true;

	/**  the maximum size of the antecedent of rules (optional) */
	int maxAntecedentSize = Integer.MAX_VALUE;

	/** the maximum size of the consequent of rules (optional) */
	int maxConsequentSize = Integer.MAX_VALUE;

	/** the decay factor used for penalty of support */
	//@Todo: Decay factor verändern + testen
	double decayFactor;

	/**
	 * Default constructor
	 */
	public AlgoRULEGROWTH() {
	}


	/**
	 * The main method to run the algorithm
	 * @param minSupport : the minimum support (percentage as a double value)
	 * @param minConfidence : the minimum confidence threshold
	 * @param input : an input file path of a sequence database
	 * @param output : a file path for writing the output file containing the seq. rules.
	 * @exception IOException if error reading/writing files
	 */
	public void runAlgorithm(double minSupport, double minConfidence, String input, String output) throws IOException {
		System.out.println("Got the following path: " + input);
		try {
			// read the input database
			database = new SequenceDatabase();
			database.loadFile(input);
		} catch (Exception e) {
			e.printStackTrace();
		}
		// convert minimum support to an absolute minimum support (double)
		this.minsupAbsolute = minSupport * database.size();
		this.minsupRelative = minSupport;
		System.out.println("Minsup absolute " + this.minsupAbsolute);
		// run the algorithm  with the just calculated absolute minimum support
		runAlgorithm(input, output, minsupAbsolute, minConfidence);
	}

	/**
	 * The main method to run the algorithm
	 * @param absoluteMinsup : the minimum support as a double value (an absolute minimum support)
	 * @param minConfidence : the minimum confidence threshold
	 * @param input : an input file path of a sequence database
	 * @param output : a file path for writing the output file containing the seq. rules.
	 * @exception IOException if error reading/writing files
	 */
	public void runAlgorithm(String input, String output, double absoluteMinsup, double minConfidence) throws IOException {
		// save the minimum confidence parameter
		this.minConfidence = minConfidence;
		// reinitialize the number of rules found
		ruleCount = 0;

		// if the database was not loaded, then load it.
		if(database == null){
			try {
				database = new SequenceDatabase();
				database.loadFile(input);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		//System.out.println(database.toString());

		// reset the stats for memory usage
		MemoryLogger.getInstance().reset();

		// prepare the object for writing the output file
		writer = new BufferedWriter(new FileWriter(output));

		// if minsup is 0, set it to 1 to avoid generating
		// rules not in the database
		this.minsupAbsolute =  absoluteMinsup;
		System.out.println("Minsup absolute " + this.minsupAbsolute);
		if(this.minsupAbsolute == 0){ // protection
			this.minsupAbsolute = 1;
		}
		System.out.println("=== START ===\nMinsup: " + minsupAbsolute + "\nMinconf: " + minConfidence + "\nInput: " + input + "\nOutput: " + output);

		// save the start time
		timeStart = System.currentTimeMillis(); // for stats

		// Remove infrequent items from the database in one database scan.
		// Then perform another database scan to count the
		// the support of each item in the same database scan
		// and their occurrences.
		removeItemsThatAreNotFrequent(database);

		// Put frequent items in a list.
		List<Integer> listFrequents = new ArrayList<Integer>();
		// for each item
		for(Entry<Integer,Map<Integer, Occurence>> entry : mapItemCount.entrySet()){
			// if it is frequent
			if(entry.getValue().size() >= minsupAbsolute){
				// add it to the list
				listFrequents.add(entry.getKey());
			}
		}
		// ======================================== fine

		// We will now try to generate rules with one item in the
		// antecedent and one item in the consequent using
		// the frequent items.
		// For each pair of frequent items i and j such that i != j
		for(int i=0; i< listFrequents.size(); i++){
			double mean_gapJI = .0;
			double mean_gapIJ = .0;
			// get the item I and its map of occurences
			Integer intI = listFrequents.get(i);
			Map<Integer, Occurence> occurencesI = mapItemCount.get(intI);
			// get the tidset of item I
			Set<Integer> tidsI = occurencesI.keySet();

			for(int j=i+1; j< listFrequents.size(); j++){
				// get the item j and its map of occurences
				Integer intJ = listFrequents.get(j);
				//System.out.println("Checking if " + intI + " and " + intJ + " are frequent together");
				Map<Integer,Occurence> occurencesJ = mapItemCount.get(intJ);
				// get the tidset of item J
				Set<Integer> tidsJ = occurencesJ.keySet();

				// (1) We will now calculate the tidsets
				// of I-->J  and the rule J-->I.

				// initialize the sets
				Set<Integer> tidsIJ = new HashSet<Integer>(); // tidset of I-->J
				Set<Integer> tidsJI = new HashSet<Integer>(); // tidset of J-->I

				// save all minimum gap distances of IJ in a list (JI respectively)
				ArrayList<Double> distsIJ = new ArrayList<>();
				ArrayList<Double> distsJI = new ArrayList<>();

				// for each occurence of I
				for(Entry<Integer, Occurence> entryOccI : occurencesI.entrySet()){
					// get the occurence of J in the same sequence
					Occurence occJ = occurencesJ.get(entryOccI.getKey());
					//System.out.println("Therefore checking tid " + entryOccI.getKey());
					// if J appears in that sequence
					if(occJ !=  null){
						ArrayList<Short> occsJ = occJ.allItemsets;
						ArrayList<Short> occsI = entryOccI.getValue().allItemsets;

						ArrayList<Double> temp = new ArrayList<>();
						// if J appeared before I in that sequence,
						// then we put this tid in the tidset of  J-->I
						if(occJ.firstItemset < entryOccI.getValue().lastItemset){
							//System.out.println("Adding to the list of J-->I because firstItemset of J appears before last Itemset of I");
							tidsJI.add(entryOccI.getKey());
							for(Short itemsetIndexJ : occsJ){
								for (Short itemsetIndexI : occsI){
									//System.out.println("ItemsetIndexI: " + itemsetIndexI);
									//System.out.println("ItemsetIndexJ: " + itemsetIndexJ);
									if(itemsetIndexJ < itemsetIndexI){
										temp.add(Math.pow(decayFactor, (itemsetIndexI-itemsetIndexJ)-1));
										//System.out.println("Found gap of " + ((itemsetIndexI-itemsetIndexJ)-1));
									}
								}
							}
						}
						if(!temp.isEmpty()){
							double min_gap = Collections.max(temp);
							distsJI.add(min_gap); // take decay resulting from least gapsize
							mean_gapJI += min_gap;
						}

						temp = new ArrayList<>();
						// if I appeared before J in that sequence,
						// then we put this tid in the tidset of  I-->J
						if(entryOccI.getValue().firstItemset < occJ.lastItemset){
							//System.out.println("Adding to the list of I-->J because firstItemset of I appears before last Itemset of J");
							tidsIJ.add(entryOccI.getKey());
							for(Short itemsetIndexI : occsI){
								for (Short itemsetIndexJ : occsJ){
									//System.out.println("ItemsetIndexI: " + itemsetIndexI);
									//System.out.println("ItemsetIndexJ: " + itemsetIndexJ);
									if(itemsetIndexI < itemsetIndexJ){
										temp.add(Math.pow(decayFactor, (itemsetIndexJ-itemsetIndexI)-1));
										//System.out.println("Found gap of " + ((itemsetIndexJ-itemsetIndexI)-1));
									}
								}
							}
						}
						if(!temp.isEmpty()){
							double min_gap = Collections.max(temp);
							distsIJ.add(min_gap); // take decay resulting from least gapsize
							mean_gapIJ += min_gap;
						}
					}
				}

				float minDecayedSupIJ = 0;
				for(Double d : distsIJ) {
					minDecayedSupIJ += d;
				}
				if(distsIJ.size() > 0){
					mean_gapIJ /= distsIJ.size();
				} else {
					mean_gapIJ = -1;

				}
				//minDecayedSupIJ /= distsIJ.size();

				float minDecayedSupJI = 0;
				for(Double d : distsJI) {
					minDecayedSupJI += d;
				}
				if(distsIJ.size() > 0){
					mean_gapJI /= distsJI.size();
				} else {
					mean_gapJI = -1;

				}
				//minDecayedSupJI /= distsJI.size();

//				System.out.println("################ Checking rule: [" + intI + "]-->[" + intJ + "]");
//				System.out.println("MinSupp: " + minsupAbsolute);
//				System.out.println("SuppIJ: " + tidsIJ.size());
//				System.out.println("SuppJI: " + tidsJI.size());
//				System.out.println("DecaySuppIJ: " + minDecayedSupIJ);
//				System.out.println("DecaySuppJI: " + minDecayedSupJI);

				// (2) check if the two itemsets have enough common tids
				// if not, we don't need to generate a rule for them.

				// create itemset of the rule I ==> J
				int[] itemsetI = new int[1];
				itemsetI[0]= intI;
				int[] itemsetJ = new int[1];
				itemsetJ[0]= intJ;

				// create rule IJ
				//if(tidsIJ.size() >= minsupAbsolute){
//				System.out.println("DecaySuppIJ >= minsupAbsolute: " + minDecayedSupIJ + ">=" + minsupAbsolute);
				if(minDecayedSupIJ >= minsupAbsolute){
					// calculate the confidence of I ==> J
					double confIJ = ((double)tidsIJ.size()) / occurencesI.size();

//					System.out.println("Confidence of [" + intI + "]-->[" + intJ + "]: " + confIJ);



					// if the confidence is high enough, save the rule
					if(confIJ >= minConfidence){
						saveRule(minDecayedSupIJ, tidsIJ, confIJ, itemsetI, itemsetJ, mean_gapIJ);
						if(DEBUG) {
							Rule rule = new Rule(itemsetI, itemsetJ, tidsI, tidsJ, tidsIJ, occurencesI, occurencesJ);
							allRulesFoundForDEBUG.add(rule);
						}
					}
				}
				// recursive call to try to expand the rule on the left and
				// right sides
				// outside the if condition because anti-monotonicity can not be applied in this case
				// although pruning strategies exist, where anti-monotonicity applies --> procedure has to be changes
				if(itemsetI.length < maxAntecedentSize) {
					//System.out.println("Expanding left with I: " + Arrays.toString(itemsetI) + " and J: " + Arrays.toString(itemsetJ));
					expandLeft(itemsetI, itemsetJ, tidsI, tidsIJ, occurencesJ);
				}
				if(itemsetJ.length < maxConsequentSize) {
					//System.out.println("Expanding right with I: " + Arrays.toString(itemsetI) + " and J: " + Arrays.toString(itemsetJ));
					expandRight(itemsetI, itemsetJ, tidsI, tidsJ, tidsIJ, occurencesI, occurencesJ);
				}


				// check if J ==> I has enough common tids
				// If yes, we create the rule J ==> I
				//if(tidsJI.size() >= minsupAbsolute){

				// create itemset of the rule J ==> I
				itemsetI = new int[1];
				itemsetI[0]= intI;
				itemsetJ = new int[1];
				itemsetJ[0]= intJ;

//				System.out.println("DecaySuppJI >= minsupAbsolute: " + minDecayedSupJI + ">=" + minsupAbsolute);
				if(minDecayedSupJI >= minsupAbsolute){
					// calculate the confidence
					double confJI = ((double)tidsJI.size()) / occurencesJ.size();

//					System.out.println("Confidence of [" + intJ + "]-->[" + intI + "]: " + confJI);

					// if the confidence is high enough, save the rule
					if(confJI >= minConfidence){
						saveRule(minDecayedSupJI, tidsJI, confJI, itemsetJ, itemsetI, mean_gapJI);
						if(DEBUG) {
							Rule rule = new Rule(itemsetJ, itemsetI, tidsJ,  tidsI, tidsJI, occurencesJ, occurencesI);
							allRulesFoundForDEBUG.add(rule);
						}
					}

				}
				// recursive call to try to expand the rule on the left and
				// right sides
				if(itemsetI.length < maxConsequentSize) {
					expandRight(itemsetJ, itemsetI, tidsJ,  tidsI, tidsJI, occurencesJ, occurencesI);
				}
				if(itemsetJ.length < maxAntecedentSize) {
					expandLeft(itemsetJ, itemsetI, tidsJ, tidsJI, occurencesI);
				}
			}
		}

		// CHECK FOR REDUNDANT RULES
		if(DEBUG) {
			for(int i=0; i < allRulesFoundForDEBUG.size(); i++) {
				for(int j=i+1; j < allRulesFoundForDEBUG.size(); j++) {
					Rule rule1 = allRulesFoundForDEBUG.get(i);
					Rule rule2 = allRulesFoundForDEBUG.get(j);
					Arrays.sort(rule1.itemsetI);
					Arrays.sort(rule1.itemsetJ);
					Arrays.sort(rule2.itemsetI);
					Arrays.sort(rule2.itemsetJ);
					if(Arrays.equals(rule1.itemsetI, rule2.itemsetI) &&
							Arrays.equals(rule1.itemsetJ, rule2.itemsetJ)) {
						throw new RuntimeException(" DUPLICATE RULES FOUND");
					}
				}
			}

		}
		// save end time
		timeEnd = System.currentTimeMillis();

		// close the file
		writer.close();

		// after the algorithm ends, we don't need a reference to the database anymore.
		database = null;
		printStats();
	}

	/**
	 * Save a rule I ==> J to the output file
	 * @param tidsIJ the tids containing the rule
	 * @param confIJ the confidence
	 * @param itemsetI the left part of the rule
	 * @param itemsetJ the right part of the rule
	 * @throws IOException exception if error writing the file
	 */
	private void saveRule(double decayedSupIJ, Set<Integer> tidsIJ, double confIJ, int[] itemsetI, int[] itemsetJ, double mean_gap) throws IOException {
		//System.out.println("######### SAVING NEW RULE #########");

		// increase the number of rule found
		ruleCount++;
//
		Arrays.sort(itemsetI);
		Arrays.sort(itemsetJ);
		//System.out.println(Arrays.toString(itemsetI) + " ==> " + Arrays.toString(itemsetJ) + " sup: " + tidsIJ.size() + "  conf : " + confIJ);
//
		// create a string buffer
		StringBuilder buffer = new StringBuilder();

		// write itemset 1 (antecedent)
		for(int i=0; i<itemsetI.length; i++){
			buffer.append(itemsetI[i]);
			if(i != itemsetI.length -1){
				buffer.append(",");
			}
		}

		// write separator
		buffer.append(" ==> ");

		// write itemset 2  (consequent)
		for(int i=0; i<itemsetJ.length; i++){
			buffer.append(itemsetJ[i]);
			if(i != itemsetJ.length -1){
				buffer.append(",");
			}
		}
		// write support
		buffer.append(" #SUP: ");
		buffer.append(tidsIJ.size());
		// write decayed support
		buffer.append(" #DECSUP: ");
		buffer.append(decayedSupIJ);
		// write confidence
		buffer.append(" #CONF: ");
		buffer.append(confIJ);
		// write mean_gap
		buffer.append(" #MEANGAP: ");
		buffer.append(mean_gap);

		//System.out.println("\n\t###############################");
		//System.out.println("\t########## SAVING RULE ########");
		//System.out.println("\t###############################\n");
		writer.write(buffer.toString());
//		System.out.println(buffer.toString());
		writer.newLine();
	}


	/**
	 * This method search for items for expanding left side of a rule I --> J
	 * with any item c. This results in rules of the form I U{c} --> J. The method makes sure that:
	 *   - c  is not already included in I or J
	 *   - c appear at least minsup time in tidsIJ before last occurence of J
	 *   - c is lexically bigger than all items in I
	 * @throws IOException
	 */
	private void expandLeft(int [] itemsetI, int[] itemsetJ, Collection<Integer> tidsI,
							Collection<Integer> tidsIJ,
							Map<Integer, Occurence> occurencesJ) throws IOException {
//		System.out.println("EXPANDLEFT");
//		System.out.println("Trying to expand " + Arrays.toString(itemsetI) + " ==> " + Arrays.toString(itemsetJ));
		// The following map will be used to count the support of each item
		// c that could potentially extend the rule.
		// The map associated a set of tids (value) to an item (key).
		Map<Integer, Set<Integer>> frequentItemsC  = new HashMap<Integer, Set<Integer>>();

		// We scan the sequence where I-->J appear to search for items c
		// that we could add to generate a larger rule  IU{c} --> J
		int left = tidsIJ.size();  // the number of tid containing I-->J

		// For each tid of sequence containing I-->J
		for(Integer tid : tidsIJ){
			// get the sequence and occurences of J in that sequence
			Sequence sequence = database.getSequences().get(tid);
			Occurence end = occurencesJ.get(tid);

			// for each itemset before the last occurence of J in that sequence
			itemLoop:	for(int k=0; k < end.lastItemset; k++){
				List<Integer> itemset = sequence.get(k);
				// for each item c in that itemset
				for(int m=0; m< itemset.size(); m++){
					Integer itemC = itemset.get(m);

					// We will consider if we could create a rule IU{c} --> J
					// If lexical order is not respected or c is included in the rule already,
					// then we cannot so return.
					if(containsLEXPlus(itemsetI, itemC) ||  containsLEX(itemsetJ, itemC)){
						continue;
					}

					// Otherwise, we get the tidset of "c"
					Set<Integer> tidsItemC = frequentItemsC.get(itemC);

					// if this set is null, which means that "c" was not seen yet
					// when scanning the sequences from I==>J

					// these conditions also apply in case of decaySupport
					// because the best case is a gapsize = 0 which means the
					// support of a rule in a sequence results in 1 (decayFactor^{gapsize})
					if(tidsItemC == null){
						// if there is less tids left in the tidset of I-->J to be scanned than
						// the minsup, we don't consider c anymore because  IU{c} --> J
						// could not be frequent
						if(left < minsupAbsolute){
							continue itemLoop;
						}
						// if "c" was seen before but there is not enough sequences left to be scanned
						// to allow IU{c} --> J to reach the minimum support threshold
					}else if(tidsItemC.size() + left < minsupAbsolute){
						// remove c and continue the loop of items
						tidsItemC.remove(itemC);
						continue itemLoop;
					}

					// otherwise, if we did not see "c" yet, create a new tidset for "c"
					if(tidsItemC == null){
						tidsItemC = new HashSet<Integer>(tidsIJ.size());
						frequentItemsC.put(itemC, tidsItemC);
					}
					// add the current tid to the tidset of "c"
					tidsItemC.add(tid);
				}
			}
			left--;  // decrease the number of sequences left to be scanned
		}

//		System.out.println("Found the following frequent items with which expansion possible: " + frequentItemsC.toString());
		// For each item c found, we create a rule	IU{c} ==> J
		for(Entry<Integer, Set<Integer>> entry : frequentItemsC.entrySet()){
			double mean_gap = .0;
			Integer itemC = entry.getKey(); // item
//			System.out.println("Using following item to expand: " + itemC);
			// get the tidset IU{c} ==> J
			Set<Integer> tidsIC_J = entry.getValue(); // tids in which IUc-->J occurs

			//System.out.println(tidsIC_J.size() + " transaction(s) has/have the expanded rule, namely: " + tidsIC_J.toString());

			ArrayList<Double> minDistsIC_J = new ArrayList<>();
			for(Integer tid : tidsIC_J) {
				int minDistIC_J = Integer.MAX_VALUE;
				//System.out.println("Looking at TID: " + tid + "\nand check whether expansion itemsetI: " +
				//		Arrays.toString(itemsetI) + "\nwith " + itemC + " works out.\nItemsetJ: " + Arrays.toString(itemsetJ) +
				//		".\nSequence: " + database.getSequences().get(tid));
				ArrayList<Short> occsI = new ArrayList<>();
				ArrayList<Short> occsJ = new ArrayList<>();
				for (int item : itemsetI) {
					Occurence o = mapItemCount.get(item).get(tid);
					occsI.addAll(o.allItemsets);
				}
				ArrayList<Short> occsC = new ArrayList<>(mapItemCount.get(itemC).get(tid).allItemsets);

				for (int item : itemsetJ) {
					occsJ.addAll(mapItemCount.get(item).get(tid).allItemsets);
				}
				Collections.sort(occsI);
				Collections.sort(occsC);
				Collections.sort(occsJ);
				//System.out.println("Occs I: " + occsI);
				//System.out.println("Occs C: " + occsC);
				//System.out.println("Occs J: " + occsJ);

				for(Short occJ : occsJ){
					int distIC_J = Integer.MAX_VALUE;
					// Check minimum distance between I and J
					int minDistI = Integer.MAX_VALUE;
					for(Short x : occsI){
						int temp = occJ-x;
						if (temp >= 0 && temp < minDistI){
							minDistI = temp;
						}
					}
					// Check minimum distance between C and J
					int minDistC = Integer.MAX_VALUE;
					for(Short x : occsC){
						int temp = occJ-x;
						if (temp >= 0 && temp < minDistC){
							minDistC = temp;
						}
					}
					// If I and C both occur before J take the minimum distance of both
					if (minDistC < Integer.MAX_VALUE && minDistI < Integer.MAX_VALUE){
						if (minDistC < minDistI) {
							distIC_J = minDistC;
						} else {
							distIC_J = minDistI;
						}
					}
					// Search for minimum distance between IUc and J over the whole sequence
					if (distIC_J < minDistIC_J){
						minDistIC_J = distIC_J;
					}
				}
				// Since I and c both have to occur at least once before J we have the minimum gapsize
				// over all the sequences where IUc-->J occurs
				minDistsIC_J.add(Math.pow(decayFactor, minDistIC_J));
				mean_gap += minDistIC_J;
			}

			double decayedSupp = .0;
			for(Double d : minDistsIC_J){
				decayedSupp += d;

			}
			if(minDistsIC_J.size() > 0){
				mean_gap /= minDistsIC_J.size();
			} else {
				mean_gap = -1;

			}
			//mean_gap /= minDistsIC_J.size();
			//avgDecayedSupp /= minDistsIC_J.size();

			//System.out.println("MinSupp: " + minsupAbsolute);
			//System.out.println("SuppIC_J: " + tidsIC_J.size());
			//System.out.println("DecaySuppIC_J: " + decayedSupp);


			// ========================== This has to be done independently of the check if the rule is frequent because
			// anti-monotonicity does not hold with this implementation
			// Calculate tids containing IU{c} which is necessary
			// to calculate the confidence
			Set<Integer> tidsIC = new HashSet<Integer>(tidsI.size());
			for(Integer tid: tidsI){
				if(mapItemCount.get(itemC).containsKey(tid)){
					tidsIC.add(tid);
				}
			}
			// create the itemset IU{c}
			int [] itemsetIC = new int[itemsetI.length+1];
			System.arraycopy(itemsetI, 0, itemsetIC, 0, itemsetI.length);
			itemsetIC[itemsetI.length] = itemC;

//			System.out.println("ExpandLeft: DecaySupp >= minsupAbsolute: " + decayedSupp + ">=" + minsupAbsolute);
			// if the support of IU{c} ==> J is enough
			//if(tidsIC_J.size() >= minsupAbsolute){
			if(decayedSupp >= minsupAbsolute){
				// Create rule and calculate its confidence of IU{c} ==> J
				// defined as:  sup(IU{c} -->J) /  sup(IU{c})
				double confIC_J = ((double)tidsIC_J.size()) / tidsIC.size();

				// if the confidence is high enough, then it is a valid rule
				if(confIC_J >= minConfidence){
					// save the rule
					saveRule(decayedSupp, tidsIC_J, confIC_J, itemsetIC, itemsetJ, mean_gap);
					if(DEBUG) {
						Rule newRule = new Rule(itemsetIC, itemsetJ, tidsIC, null, tidsIC_J, null, occurencesJ);
						allRulesFoundForDEBUG.add(newRule);
					}
				}
			}
			// recursive call to expand left side of the rule
			if(itemsetI.length < maxAntecedentSize) {
				expandLeft(itemsetIC, itemsetJ, tidsIC, tidsIC_J, occurencesJ);
			}
			// ===============================
		}
		// check the memory usage
		MemoryLogger.getInstance().checkMemory();
	}

	/**
	 * This method search for items for expanding left side of a rule I --> J
	 * with any item c. This results in rules of the form I --> J U�{c}. The method makes sure that:
	 *   - c  is not already included in I or J
	 *   - c appear at least minsup time in tidsIJ after the first occurence of I
	 *   - c is lexically bigger than all items in J
	 * @throws IOException
	 */
	private void expandRight(int[] itemsetI, int[] itemsetJ,
							 Set<Integer> tidsI,
							 Collection<Integer> tidsJ,
							 Collection<Integer> tidsIJ,
							 Map<Integer, Occurence> occurencesI,
							 Map<Integer, Occurence> occurencesJ) throws IOException {
//		if(true)
//    	return;
//		System.out.println("EXPANDRIGHT");
//		System.out.println("Trying to expand " + Arrays.toString(itemsetI) + " ==> " + Arrays.toString(itemsetJ));
		// The following map will be used to count the support of each item
		// c that could potentially extend the rule.
		// The map associated a set of tids (value) to an item (key).
		Map<Integer, Set<Integer>> frequentItemsC  = new HashMap<Integer, Set<Integer>>();

		// we scan the sequence where I-->J appear to search for items c that we could add.
		// for each sequence containing I-->J.
		int left = tidsIJ.size();

		// For each tid of sequence containing I-->J
		for(Integer tid : tidsIJ){
			// get the sequence and get occurences of I in that sequence
			Sequence sequence = database.getSequences().get(tid);
//			System.out.println("Checking sequence: " + tid + " = " + sequence.toString());
			Occurence first = occurencesI.get(tid);

			// for each itemset after the first occurence of I in that sequence
			for(int k=first.firstItemset+1; k < sequence.size(); k++){
				List<Integer> itemset = sequence.get(k);
				// for each item
				itemLoop:	for(int m=0; m< itemset.size(); m++){
					// for each item c in that itemset
					Integer itemC = itemset.get(m);

					// We will consider if we could create a rule I --> J U{c}
					// If lexical order is not respected or c is included in the rule already,
					// then we cannot so the algorithm return.
					if(containsLEX(itemsetI, itemC) ||  containsLEXPlus(itemsetJ, itemC)){
						continue;
					}
					Set<Integer> tidsItemC = frequentItemsC.get(itemC);
					System.out.println("" + itemC + " " + tidsItemC + " " + left);
					// if "c" was seen before but there is not enough sequences left to be scanned
					// to allow IU --> J {c} to reach the minimum support threshold

					// these conditions also apply in case of decaySupport
					// because the best case is a gapsize = 0 which means the
					// support of a rule in a sequence results in 1 (decayFactor^{gapsize})

					if(tidsItemC == null){
						if(left < minsupAbsolute){
							continue itemLoop;
						}
					}else if(tidsItemC.size() + left < minsupAbsolute){
						// if "c" was seen before but there is not enough sequences left to be scanned
						// to allow I--> JU{c}  to reach the minimum support threshold,
						// remove "c" and continue the loop of items
						tidsItemC.remove(itemC);
						continue itemLoop;
					}

					if(tidsItemC == null){
						// otherwise, if we did not see "c" yet, create a new tidset for "c"
						tidsItemC = new HashSet<Integer>(tidsIJ.size());
						frequentItemsC.put(itemC, tidsItemC);
					}
					// add the current tid to the tidset of "c"
					tidsItemC.add(tid);
				}
			}
			left--;  // decrease the number of sequences left to be scanned
		}
//		System.out.println("Found the following frequent items with which expansion possible: " + frequentItemsC.toString());
		// For each item c found, we create a rule	I ==> JU {c}
		for(Entry<Integer, Set<Integer>> entry : frequentItemsC.entrySet()){
			double mean_gap = .0;
			Integer itemC = entry.getKey();
//			System.out.println("Using following item to expand: " + itemC);
			// get the tidset of I ==> JU {c}
			Set<Integer> tidsI_JC = entry.getValue();
			//System.out.println(tidsI_JC.size() + " transaction(s) has/have the expanded rule, namely: " + tidsI_JC.toString());

			ArrayList<Double> minDistsI_JC = new ArrayList<>();
			for(Integer tid : tidsI_JC) {
//				System.out.println("Current sequence " + tid);
				int minDistI_JC = Integer.MAX_VALUE;
				//System.out.println("Looking at TID: " + tid + "\nand check whether expansion itemsetJ: " +
				//		Arrays.toString(itemsetJ) + "\nwith " + itemC + " works out.\nItemsetJ: " + Arrays.toString(itemsetI) +
				//		".\nSequence: " + database.getSequences().get(tid));
				ArrayList<Short> occsI = new ArrayList<>();
				ArrayList<Short> occsJ = new ArrayList<>();
				for (int item : itemsetI) {
					Occurence o = mapItemCount.get(item).get(tid);
					occsI.addAll(o.allItemsets);
				}
				ArrayList<Short> occsC = new ArrayList<>(mapItemCount.get(itemC).get(tid).allItemsets);
//				System.out.println("ItemsetJ: " + Arrays.toString(itemsetJ));
				for (int item : itemsetJ) {
//					System.out.println("Get mapping of tid to occurence: " + mapItemCount.get(item).get(tid));
//					System.out.println("Get mapping of item to tid and occurence: " + mapItemCount.get(item));
//					System.out.println("Get all occurences of item in sequence: " + mapItemCount.get(item).get(tid).allItemsets);
					occsJ.addAll(mapItemCount.get(item).get(tid).allItemsets);
				}
				Collections.sort(occsI);
				Collections.sort(occsC);
				Collections.sort(occsJ);
//				System.out.println("Occs I: " + occsI);
//				System.out.println("Occs J: " + occsJ);
//				System.out.println("Occs C: " + occsC);

				for(Short occI : occsI){
					int distI_JC = Integer.MAX_VALUE;
					// Check minimum distance between I and J
					int minDistJ = Integer.MAX_VALUE;
					for(Short x : occsJ){
						int temp = x-occI;
						if (temp >= 0 && temp < minDistJ){
							minDistJ = temp;
						}
					}
					// Check minimum distance between C and J
					int minDistC = Integer.MAX_VALUE;
					for(Short x : occsC){
						int temp = x-occI;
						if (temp >= 0 && temp < minDistC){
							minDistC = temp;
						}
					}
					// If I and C both occur before J take the minimum distance of both
					if (minDistC < Integer.MAX_VALUE && minDistJ < Integer.MAX_VALUE){
						if (minDistC < minDistJ) {
							distI_JC = minDistC;
						} else {
							distI_JC = minDistJ;
						}
					}
					// Search for minimum distance between I and JUc over the whole sequence
					if (distI_JC < minDistI_JC){
						minDistI_JC = distI_JC;
					}
				}
				// Since I and c both have to occur at least once before J we have the minimum gapsize
				// over all the sequences where IUc-->J occurs
				minDistsI_JC.add(Math.pow(decayFactor, minDistI_JC));
				mean_gap += minDistI_JC;
			}


			double decayedSupp = .0;
			for(Double d : minDistsI_JC)
				decayedSupp += d;
			//decayedSupp /= minDistsI_JC.size();

			if(minDistsI_JC.size() > 0){
				mean_gap /= minDistsI_JC.size();
			} else {
				mean_gap = -1;

			}

			//mean_gap /= minDistsI_JC.size();
//			System.out.println("MinSupp: " + minsupAbsolute);
//			System.out.println("SuppI_JC: " + tidsI_JC.size());
//			System.out.println("DecaySuppI_JC: " + decayedSupp);

			// ========================== This has to be done independently of the check if the rule is frequent because
			// anti-monotonicity does not hold with this implementation
			int[] itemsetJC = new int[itemsetJ.length+1];
			System.arraycopy(itemsetJ, 0, itemsetJC, 0, itemsetJ.length);
			itemsetJC[itemsetJ.length]= itemC;

			// create the itemset JU{c} and calculate the occurences of JU{c}
			Set<Integer> tidsJC = new HashSet<Integer>(tidsJ.size());
			Map<Integer, Occurence> occurencesJC = new HashMap<Integer, Occurence>();
			// for each sequence containing J
			for(Integer tid: tidsJ){
				// Get the first and last occurences of C in that sequence
				Occurence occurenceC = mapItemCount.get(itemC).get(tid);
				// if there is an occurence
				if(occurenceC != null){
					// add the tid of the sequence to the tidset of JU{c}
					tidsJC.add(tid);
					// calculate last occurence of JU{c} depending on if
					// the last occurence of J is before the last occurence
					// of c or not.
					Occurence occurenceJ = occurencesJ.get(tid);
					if(occurenceC.lastItemset < occurenceJ.lastItemset){
						occurencesJC.put(tid, occurenceC);
					}else{
						occurencesJC.put(tid, occurenceJ);
					}
				}
			}

//			System.out.println("ExpandRight: DecaySupp >= minsupAbsolute: " + decayedSupp + ">=" + minsupAbsolute);
			// if the support of I ==> JU{c} is enough
			//if(tidsI_JC.size() >= minsupAbsolute){
			if(decayedSupp >= minsupAbsolute){

				// Create rule I ==> J U{c} and calculate its confidence
				// defined as:  sup(I -->J U{c}) /  sup(I)
				double confI_JC = ((double)tidsI_JC.size()) / tidsI.size();

				// if the confidence is enough
				if(confI_JC >= minConfidence){
					// then it is a valid rule so save it
					saveRule(decayedSupp, tidsI_JC, confI_JC, itemsetI, itemsetJC, mean_gap);
					if(DEBUG) {
						Rule newRule = new Rule(itemsetI, itemsetJC, tidsI, tidsJC, tidsI_JC, occurencesI, occurencesJC);
						allRulesFoundForDEBUG.add(newRule);
					}
				}
			}
			// recursively try to expand the left and right side
			// of the rule
			if(itemsetJC.length < maxConsequentSize) {
				expandRight(itemsetI, itemsetJC, tidsI, tidsJC, tidsI_JC, occurencesI, occurencesJC);  // occurencesJ
			}
			if(itemsetI.length < maxAntecedentSize) {
				expandLeft(itemsetI, itemsetJC,  tidsI, tidsI_JC, occurencesJC);  // occurencesJ
			}
			//========================
		}
		// check the memory usage
		MemoryLogger.getInstance().checkMemory();
	}


	/**
	 * This method calculate the frequency of each item in one database pass.
	 * Then it remove all items that are not frequent.
	 * @param database : a sequence database
	 * @return A map such that key = item
	 *                         value = a map  where a key = tid  and a value = Occurence
	 * This map allows knowing the frequency of each item and their first and last occurence in each sequence.
	 */
	private Map<Integer, Map<Integer, Occurence>> removeItemsThatAreNotFrequent(SequenceDatabase database) {
		// (1) Count the support of each item in the database in one database pass
		mapItemCount = new HashMap<Integer, Map<Integer, Occurence>>(); // <item, Map<tid, occurence>>

		// for each sequence in the database
		for(int k=0; k< database.size(); k++){
			Sequence sequence = database.getSequences().get(k);
			// for each itemset in that sequence
			for(short j=0; j< sequence.getItemsets().size(); j++){
				List<Integer> itemset = sequence.get(j);
				// for each item in that itemset
				for(int i=0; i< itemset.size(); i++){
					Integer itemI = itemset.get(i);

					// get the map of occurences of that item
					Map<Integer, Occurence> occurences = mapItemCount.get(itemI);
					// if this map is null, create a new one
					if(occurences == null){
						occurences = new HashMap<Integer, Occurence>();
						mapItemCount.put(itemI, occurences);
					}
					// then update the occurence by adding j as the
					// last occurence in sequence k
					Occurence occurence = occurences.get(k);
					if(occurence == null){
						occurence = new Occurence(j, j);
						occurence.allItemsets.add(j);
						occurences.put(k, occurence);
					}else{
						occurence.lastItemset = j;
						occurence.allItemsets.add(j);
					}
				}
			}
		}
//		System.out.println("NUMBER OF DIFFERENT ITEMS : " + mapItemCount.size());
		// (2) remove all items that are not frequent from the database

		// for each sequence
		for(Sequence sequence : database.getSequences()){
			int i=0;

			// for each itemset
			while(i < sequence.getItemsets().size()){
				List<Integer> itemset = sequence.getItemsets().get(i);
				int j=0;

				// for each item
				while(j < itemset.size()){
					// if the item is not frequent remove it
					if( mapItemCount.get(itemset.get(j)).size() < minsupAbsolute){
						itemset.remove(j);
					}else{
						// otherwise go to next item
						j++;
					}
				}
				i++;  // go to next itemset
			}
		}
		// return the map of occurences of items
		return mapItemCount;
	}

	/**

	 * This method checks if the item "item" is in the itemset.
	 * It asumes that items in the itemset are sorted in lexical order
	 * This version also checks that if the item "item" was added it would be the largest one
	 * according to the lexical order.
	 * @param itemset an itemset
	 * @param item  the item
	 * @return return true if the above conditions are met, otherwise false
	 */
	boolean containsLEXPlus(int[] itemset, int item) {
		// for each item in itemset
		for(int i=0; i< itemset.length; i++){
			// check if the current item is equal to the one that is searched
			if(itemset[i] == item){
				// if yes return true
				return true;
				// if the current item is larger than the item that is searched,
				// then return true because if if the item "item" was added it would be the largest one
				// according to the lexical order.
			}else if(itemset[i] > item){
				return true; // <-- XXXX
			}
		}
		// if the searched item was not found, return false.
		return false;
	}

	/**
	 * This method checks if the item "item" is in the itemset.
	 * It assumes that items in the itemset are sorted in lexical order
	 * @param itemset an itemset
	 * @param item  the item
	 * @return return true if the item
	 */
	boolean containsLEX(int[] itemset, int item) {
		// for each item in itemset
		for(int i=0; i< itemset.length; i++){
			// check if the current item is equal to the one that is searched
			if(itemset[i] == item){
				// if yes return true
				return true;
				// if the current item is larger than the item that is searched,
				// then return false because of the lexical order.
			}else if(itemset[i] > item){
				return false;  // <-- xxxx
			}
		}
		// if the searched item was not found, return false.
		return false;
	}

	/**
	 * Set the number of items that a rule antecedent should contain (optional).
	 * @param maxAntecedentSize the maximum number of items
	 */
	public void setMaxAntecedentSize(int maxAntecedentSize) {
		this.maxAntecedentSize = maxAntecedentSize;
	}

	public void setDecayFactor(double decayFactor) { this.decayFactor = decayFactor;}
	/**
	 * Set the number of items that a rule consequent should contain (optional).
	 * @param maxConsequentSize the maximum number of items
	 */
	public void setMaxConsequentSize(int maxConsequentSize) {
		this.maxConsequentSize = maxConsequentSize;
	}

	/**
	 * Print statistics about the last algorithm execution to System.out.
	 */
	public void printStats() {
		System.out.println("=============  RULEGROWTH - STATS ========");
		System.out.println("Sequential rules count: " + ruleCount);
		System.out.println("Total time: " + (timeEnd - timeStart) + " ms");
		System.out.println("Max memory: " + MemoryLogger.getInstance().getMaxMemory());
		System.out.println("==========================================");
	}

	public static void main (String[] args) {
		AlgoRULEGROWTH rg = new AlgoRULEGROWTH();
		//read in all arguments:
		double sup = Double.parseDouble(args[0]);
		double conf = Double.parseDouble(args[1]);
		String input = args[2];
		String output = args[3];
		rg.setMaxConsequentSize(Integer.parseInt(args[4]));
		rg.setDecayFactor(Double.parseDouble(args[5]));
//		System.out.println("Received the following input");
//		for (String s  : args){
//			System.out.println(s);
//		}
		System.out.println("Sup " + sup + " Conf " + conf);
		try {
			rg.runAlgorithm(sup, conf, input, output);
			//rg.runAlgorithm(input, output, sup, conf);
		} catch (Exception e) {
			System.out.println("ERROR in main: " + e);
		}
	}
}