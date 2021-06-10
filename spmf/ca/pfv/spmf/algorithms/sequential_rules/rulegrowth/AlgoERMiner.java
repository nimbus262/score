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
 * This is the original implementation of the ERMiner algorithm for mining sequential rules
 * common to several sequences where antecedent and consequent are unordered itemsets. 
 * <br/><br/>
 * The main method of this algorithm is "runAlgorithm". It output the result to a file.
 *
 * @see Occurence
 * @see Sequence
 * @see SequenceDatabase
 * @author Philippe Fournier-Viger
 */
public class AlgoERMiner {
	//*** for statistics ***/
	/** start time of latest execution */
	long timeStart = 0;  
	
	/** end time of latest execution */
	long timeEnd = 0;  
	
	/**  number of rules generated */
	int ruleCount; 
	
	//*** parameters ***/
	/** minimum confidence */
	double minConfidence;
	
	/** minimum support */
	double minsupAbsolute;

	/** minimum relative support */
	double minsupRelative;
	
	/** this is the sequence database */
	SequenceDatabase database;
	
	//*** internal variables ***/
	/** This map contains for each item (key) a map of occurences (value).
	// The map of occurences associates to sequence ID (key), an occurence of the item (value). */
	Map<Integer,  Map<Integer, Occurence>> mapItemCount;  // item, <tid, occurence> 

	/** object to write the output file */
	BufferedWriter writer = null; 
	
	/** the expandleft store */
	ExpandLeftStore store = new ExpandLeftStore();
	
	/** the sparse matrix */
	SparseMatrix matrix = new SparseMatrix();

	/** total number of candidtaes */
	private long totalCandidateCount;
	
	/** number of pruned candidates */
	private long candidatePrunedCount;
	
	/** the maximum size of the antecedent of rules (optional) */
	
	int maxAntecedentSize = Integer.MAX_VALUE;
	/**  the maximum size of the consequent of rules (optional) */
	int maxConsequentSize = Integer.MAX_VALUE;

	/** the decay factor used for penalty of support */
	//@Todo: Decay factor verändern + testen
	double decayFactor;

	/**
	 * Default constructor
	 */
	public AlgoERMiner() {
	}


	/**
	 * The main method to run the algorithm
	 * @param minSupport : the minimum support (percentage as a double value)
	 * @param minConfidence : the minimum confidence threshold
	 * @param input : an input file path of a seq   uence database
	 * @param output : a file path for writing the output file containing the seq. rules.
	 * @exception IOException if error reading/writing files
	 */
	public void runAlgorithm(double minSupport, double minConfidence, String input, String output) throws IOException {
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
		System.out.println("=== START ===\nMinsup: " + minSupport + " (" + minsupAbsolute + ")\nDecay Factor: " + decayFactor +"\nMinconf: " + minConfidence + "\nInput: " + input + "\nOutput: " + output);

		// run the algorithm  with the just calculated absolute minimum support
		runAlgorithm(input, output, minsupAbsolute, minConfidence);
	}
	
	/**
	 * The main method to run the algorithm
	 * @param absoluteMinsup : the minimum support as an integer value (an absolute minimum support)
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
//		countPruning = 0;
		
		// if the database was not loaded, then load it.
		if(database == null){
			try {
				database = new SequenceDatabase(); 
				database.loadFile(input);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		// reset the stats for memory usage
		MemoryLogger.getInstance().reset();

		// prepare the object for writing the output file
		writer = new BufferedWriter(new FileWriter(output));
		
		// if minsup is 0, set it to 1 to avoid generating
		// rules not in the database
		this.minsupAbsolute =  absoluteMinsup;
		if(this.minsupAbsolute == 0){ // protection
			this.minsupAbsolute = 1;
		}

		// save the start time
		timeStart = System.currentTimeMillis(); // for stats

		// Remove infrequent items from the database in one database scan.
		// Then perform another database scan to count the
		// the support of each item in the same database scan 
		// and their occurrences.
		if(maxAntecedentSize > 0 && maxConsequentSize > 0){
			calculateFrequencyOfEachItem(database);
			
			// =================== CALCULATE MATRIX ==============
			// for each sequence
			generateMatrix(database);
			// =================== 
		}

		// We will now generate the equivalence classes of rules of size 1*1
		// We create two maps to store the equivalence classes.
		// The map for left-equivalence classes has a key : item on the right side, value: list of rules
		Map<Integer, LeftEquivalenceClass> mapEclassLeft = new HashMap<Integer, LeftEquivalenceClass>();
		// The map for right-equivalence classes has a key : item on the left side, value: list of rules
		Map<Integer, RightEquivalenceClass> mapEclassRight = new HashMap<Integer, RightEquivalenceClass>();
		
		// for each entry in the matrix
		// entry  =  key: item I    value: a map where key: item J and value: the support of {i, j} 
		for(Entry<Integer, Map<Integer, Integer>> entry :  matrix.getMatrix().entrySet()) {
			Integer intI = entry.getKey();

			Map<Integer, Occurence> occurencesI = mapItemCount.get(intI);
			// get the tidset of item I
			Set<Integer> tidsI = occurencesI.keySet();
			
			// for each item J that co-occurs with I
			for(Entry<Integer, Integer> entryJ : entry.getValue().entrySet()) {
				// if I and J are infrequent, we don't continue
				if(entryJ.getValue() < minsupAbsolute) {
					continue;
				}
				Integer intJ = entryJ.getKey();

				Map<Integer, Occurence> occurencesJ = mapItemCount.get(intJ);

				// (1) We will now calculate the tids of I-->J  and J-->I.
				
				// initialize the tids sets
				Set<Integer> tidsIJ = new HashSet<Integer>();  // tidset of I-->J
				Set<Integer> tidsJI = new HashSet<Integer>();  // tidset of J-->I

				// save all minimum gap distances of IJ in a list (JI respectively)
				ArrayList<Double> distsIJ = new ArrayList<>();
				ArrayList<Double> distsJI = new ArrayList<>();
				// for each occurence of I
				if(occurencesI.size() < occurencesJ.size()) {
//					System.out.println("################ Checking rule: [" + intI + "]-->[" + intJ + "]");
					calculateTidsetsIJandJI(occurencesI, occurencesJ, tidsIJ, tidsJI, distsIJ, distsJI);
				}else {
//					System.out.println("################ Checking rule: [" + intJ + "]-->[" + intI + "]");
					calculateTidsetsIJandJI(occurencesJ, occurencesI, tidsJI, tidsIJ, distsJI, distsIJ);
				}


//				System.out.println("distsIJ " + distsIJ);
				double decayedSuppIJ = 0.0;
				for(Double d : distsIJ) {
					decayedSuppIJ += d;
				}

//				System.out.println("distsJI " + distsJI);
				double decayedSuppJI = 0.0;
				for(Double d : distsJI) {
					decayedSuppJI += d;
				}

//				System.out.println("DecaySuppIJ: " + decayedSuppIJ);
//				System.out.println("DecaySuppJI: " + decayedSuppJI);
				// (2) check if the two itemsets have enough common tids
				// if not, we don't need to generate a rule for them.
				
				// create rule IJ
				//if(tidsIJ.size() >= minsupAbsolute){
				if(decayedSuppIJ >= minsupAbsolute){

					// calculate the confidence of I ==> J
					double confIJ = ((double)tidsIJ.size()) / occurencesI.size();

//					System.out.println("Confidence of [" + intI + "]-->[" + intJ + "]: " + confIJ);

					// create itemset of the rule I ==> J
					int[] itemsetI = new int[]{intI};
					int[] itemsetJ = new int[]{intJ};

					Set<Integer> tidsJ = occurencesJ.keySet();

					// if the confidence is high enough, save the rule
					if(confIJ >= minConfidence){
						saveRule(decayedSuppIJ, tidsIJ, confIJ, itemsetI, itemsetJ);
					}
					if(maxAntecedentSize > 1 || maxConsequentSize > 1){
						// register the rule in the appropriate equivalence classes
						registerRule11(intI, intJ, tidsI, tidsJ, tidsIJ, occurencesI, occurencesJ, mapEclassLeft, mapEclassRight);
					}
				}
				// check if J ==> I has enough common tids
				// If yes, we create the rule J ==> I
				//if(tidsJI.size() >= minsupAbsolute){
				if(decayedSuppJI >= minsupAbsolute){
					// create itemset of the rule J ==> I
					int[] itemsetI = new int[]{intI};
					int[] itemsetJ = new int[]{intJ};
					
					// calculate the confidence
					double confJI = ((double)tidsJI.size()) / occurencesJ.size();

//					System.out.println("Confidence of [" + intJ + "]-->[" + intI + "]: " + confJI);

					Set<Integer> tidsJ = occurencesJ.keySet();

					// if the confidence is high enough, save the rule
					if(confJI >= minConfidence){
						saveRule(decayedSuppJI, tidsJI, confJI, itemsetJ, itemsetI);
					}
					// register the rule in the appropriate equivalence classes
					if(maxAntecedentSize > 1 || maxConsequentSize > 1){
						registerRule11(intJ, intI, tidsJ,  tidsI, tidsJI, occurencesJ, occurencesI, mapEclassLeft, mapEclassRight);
					}
				}
			}
		}

		// PERFORM EXPAND LEFT FOR EACH LEFT-EQUIVALENCE CLASS OF SIZE 1-1
		if(maxAntecedentSize >1){
			for(LeftEquivalenceClass eclassLeft : mapEclassLeft.values()) {
				if(eclassLeft.rules.size() != 1) {
					Collections.sort(eclassLeft.rules, new Comparator<LeftRule>() {
						public int compare(LeftRule arg0, LeftRule arg1) {
							return arg0.itemsetI[0] - arg1.itemsetI[0];
						}});
	
					expandLeft(eclassLeft);
				}
			}
		}
		
		mapEclassLeft = null;
		
		// PERFORM EXPAND RIGHT FOR EACH RIGHT-EQUIVALENCE CLASS OF SIZE 1-1
		if(maxConsequentSize >1){
			for(RightEquivalenceClass eclassRight : mapEclassRight.values()) {
				if(eclassRight.rules.size() != 1) {
					Collections.sort(eclassRight.rules, new Comparator<RightRule>() {
						public int compare(RightRule arg0, RightRule arg1) {
							return arg0.itemsetJ[0] - arg1.itemsetJ[0];
						}});
						
					expandRight(eclassRight, true);
				}
			}
		}
		
		mapEclassRight = null;
		
		// PROCESS ALL EQUIVALENCE CLASSES WITH MORE THAN ONE ITEM IN LEFT PART FOR EXPAND LEFT AFTER RIGHT...
		for(Map<Integer, List<LeftEquivalenceClass>> map : store.getStore().values()) {
			for(List<LeftEquivalenceClass> eclassList :  map.values()) {
				for(LeftEquivalenceClass eclass : eclassList) {
					if(eclass.rules.size() != 1) {
						Collections.sort(eclass.rules, new Comparator<LeftRule>() {
							public int compare(LeftRule arg0, LeftRule arg1) {
								return arg0.itemsetI[arg0.itemsetI.length-1] - arg1.itemsetI[arg1.itemsetI.length-1];
							}});
					
						expandLeft(eclass);
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

	private void registerRule11(Integer intI, Integer intJ, Set<Integer> tidsI,
			Set<Integer> tidsJ, Set<Integer> tidsIJ,
			Map<Integer, Occurence> occurencesI,
			Map<Integer, Occurence> occurencesJ,
			Map<Integer, LeftEquivalenceClass> mapEclassLeft,
			Map<Integer, RightEquivalenceClass> mapEclassRight) {
		
		// add the rule to the left equivalence class
		LeftEquivalenceClass leftClass = mapEclassLeft.get(intJ);
		if(leftClass == null) {
			leftClass = new LeftEquivalenceClass(new int[] {intJ}, tidsJ, occurencesJ);
			mapEclassLeft.put(intJ, leftClass);
		}
		LeftRule ruleL = new LeftRule(new int[] {intI}, tidsI, tidsIJ);
		leftClass.rules.add(ruleL);
		 
		// add the rule to the right equivalence class
		RightEquivalenceClass rightclass = mapEclassRight.get(intI);
		if(rightclass == null) {
			rightclass = new RightEquivalenceClass(new int[] {intI}, tidsI, occurencesI);
			mapEclassRight.put(intI, rightclass);
		}
		RightRule ruleR = new RightRule(new int[] {intJ}, tidsJ, tidsIJ, occurencesJ);
		rightclass.rules.add(ruleR);
	}


	private void calculateTidsetsIJandJI(Map<Integer, Occurence> occurencesI,
											 Map<Integer, Occurence> occurencesJ, Set<Integer> tidsIJ, Set<Integer> tidsJI, ArrayList<Double> distsIJ, ArrayList<Double> distsJI) {



		for(Entry<Integer, Occurence> entryOccI : occurencesI.entrySet()){
			Integer tid = entryOccI.getKey();
			// get the occurence of J in the same sequence
			Occurence occJ = occurencesJ.get(tid);
//			System.out.println("Therefore checking tid " + entryOccI.getKey());
			// if J appears in that sequence
			if(occJ !=  null){
				ArrayList<Short> occsJ = occJ.allItemsets;
				ArrayList<Short> occsI = entryOccI.getValue().allItemsets;

				ArrayList<Double> temp = new ArrayList<>();

				Occurence occI = entryOccI.getValue();
//				System.out.println("Occurence I: " + occI.allItemsets);
//				System.out.println("Occurence J: " + occJ.allItemsets);
				// if J appeared before I in that sequence,
				// then we put this tid in the tidset of J-->I
				if(occJ.firstItemset < occI.lastItemset){
//					System.out.println("Adding to the list of J-->I because firstItemset of J appears before lastItemset of I");
					tidsJI.add(tid);
					for(Short itemsetIndexJ : occsJ){
						for (Short itemsetIndexI : occsI){
							if(itemsetIndexJ < itemsetIndexI){
								int gap = (itemsetIndexI-itemsetIndexJ)-1;
								temp.add(Math.pow(decayFactor, gap));
//								System.out.println("Found gap of " + gap);
							}
						}
					}
//					System.out.println("Distances in this sequence: " + temp.toString());
				}
				if(!temp.isEmpty()){
//					System.out.println("Adding to distsJI");
					distsJI.add(Collections.max(temp)); // take decay resulting from least gapsize
				}
				temp = new ArrayList<>();
				// if I appeared before J in that sequence,
				// then we put this tid in the tidset of I-->J
				if(occI.firstItemset < occJ.lastItemset){
//					System.out.println("Adding to the list of I-->J because firstItemset of I appears before lastItemset of J");
					tidsIJ.add(tid);
					for(Short itemsetIndexI : occsI){
						for (Short itemsetIndexJ : occsJ){
							if(itemsetIndexI < itemsetIndexJ){
								int gap = (itemsetIndexJ-itemsetIndexI)-1;
								temp.add(Math.pow(decayFactor, gap));
//								System.out.println("Found gap of " + gap);
							}
						}
					}
//					System.out.println("Distances in this sequence: " + temp.toString());
				}
				if(!temp.isEmpty()){
//					System.out.println("Adding to distsIJ");
					distsIJ.add(Collections.max(temp)); // take decay resulting from least gapsize
				}
			}
		}

//		System.out.println("MinSupp: " + minsupAbsolute);
//		System.out.println("SuppIJ: " + tidsIJ.size());
//		System.out.println("SuppJI: " + tidsJI.size());
	}
	
	public int[] concatenate(int [] itemset, int item) {
		int[] newItemset = new int[itemset.length+1];
		System.arraycopy(itemset, 0, newItemset, 0, itemset.length);
		newItemset[itemset.length] = item;
		return newItemset;
	}
	

	private void expandLeft(LeftEquivalenceClass eclass){
//		System.out.println("EXPANDLEFT");
//		System.out.println("Trying to expand equivalence class " + eclass.toString());

		for(int w=0; w < eclass.rules.size()-1; w++){  // IMPORTANT : SIZE -1 BECAUSE THE LAST ONE HAS NOTHING LEFT FOR COMPARISON
			LeftRule rule1 = eclass.rules.get(w);
			int d = rule1.itemsetI[rule1.itemsetI.length -1];
			
			LeftEquivalenceClass rulesForRecursion 
				= new LeftEquivalenceClass(eclass.itemsetJ, eclass.tidsJ,
						eclass.occurencesJ);

			// for each rule J != I
			for(int m=w+1; m < eclass.rules.size(); m++)	{
				LeftRule rule2 = eclass.rules.get(m);

				int c = rule2.itemsetI[rule2.itemsetI.length -1];
				
				if(matrix.getCount(c, d) < minsupAbsolute){
					candidatePrunedCount++;
					totalCandidateCount++;
					continue;
				}
				totalCandidateCount++;

				
				// CALCULATE TIDS I U {C}
				Set<Integer> tidsIC = new HashSet<Integer>(); 
				
				Map<Integer, Occurence> mapC = mapItemCount.get(c);
				
				// depending of the relative size of tids of I and tids of {c} we have two ways to calculate it
				
				if(rule1.tidsI.size() < mapC.size()) {
					// EARLY SKIP OPTIMIZATION: <<<<<********
					int remains = rule1.tidsI.size();
		 			// for each sequence containing I
					for(Integer tid: rule1.tidsI){
						// Get the first and last occurences of C in that sequence
						// if there is an occurence
			    		if(mapC.get(tid) != null){
			    			// add the tid of the sequence to the tidset of IU{c}
			    			tidsIC.add(tid);
			    		}
			    		remains--;
			    		// early skip
			    		if(tidsIC.size() + remains < minsupAbsolute) {
			    			break;
			    		}
			    	}
				}else {
					int remains = mapC.size();
					// for each sequence containing I
					for(Integer tid: mapC.keySet()){
						// Get the first and last occurences of C in that sequence
						// if there is an occurence
			    		if(rule1.tidsI.contains(tid)){
			    			// add the tid of the sequence to the tidset of JU{c}
			    			tidsIC.add(tid);
			    		}
			    		remains--;
			    		if(tidsIC.size() + remains < minsupAbsolute) {
			    			break;
			    		}
			    	}
				}
				
				// CALCULATE TIDS IC ==> J
				Set<Integer> tidsIC_J = new HashSet<Integer>();
		    	
				// depending on the relative size of tids of IUJ and tids of {c} we have two ways to calculate it
				if(rule1.tidsIJ.size() < mapC.size()) {
					// for each sequence containing I
					for(Integer tid: rule1.tidsIJ){
						// Get the first and last occurences of C in that sequence
						Occurence occurenceC = mapC.get(tid);
						// if there is an occurence
			    		if(occurenceC != null){
							Occurence occurenceJ = eclass.occurencesJ.get(tid);
			    			if(occurenceC.firstItemset < occurenceJ.lastItemset){
			    				// add the tid of the sequence to the tidset of JU{c}
			    				tidsIC_J.add(tid);
				    		}
			    		}
			    	}
				}else {
					// for each sequence containing I
					for(Entry<Integer, Occurence> entryC: mapC.entrySet()){
						int tid = entryC.getKey();
						if(rule1.tidsIJ.contains(tid)) {
							// Get the first and last occurences of C in that sequence
							Occurence occurenceC = entryC.getValue();
							// if there is an occurence
							Occurence occurenceJ = eclass.occurencesJ.get(tid);
			    			if(occurenceC.firstItemset < occurenceJ.lastItemset){
				    			// add the tid of the sequence to the tidset of JU{c}
				    			tidsIC_J.add(tid);
			    			}
						}
			    	}
				}

				ArrayList<Double> minDistsIC_J = getMinDists(rule1.itemsetI, eclass.itemsetJ, tidsIC_J, c,"left");
				double decayedSupp = .0;
				for(Double dist : minDistsIC_J) decayedSupp += dist;

				//if(tidsIC_J.size() >= minsupAbsolute) {
				if(decayedSupp >= minsupAbsolute) {
					// Create rule and calculate its confidence of IU{c} ==> J 
			    	// defined as:  sup(IU{c} -->J) /  sup(IU{c})			
					double confIC_J = ((double) tidsIC_J.size()) / tidsIC.size();
		
					// try to combine the rules
					int itemsetIC[] = new int[rule1.itemsetI.length+1];
					System.arraycopy(rule1.itemsetI, 0, itemsetIC, 0, rule1.itemsetI.length);
					itemsetIC[rule1.itemsetI.length] = c;

					LeftRule newRule = new LeftRule(itemsetIC, tidsIC, tidsIC_J);
					
					// if the confidence is high enough, then it is a valid rule
					if(confIC_J >= minConfidence){
						// save the rule
						saveRule(decayedSupp, tidsIC_J, confIC_J, itemsetIC,  eclass.itemsetJ);
					}
		
					if(newRule.itemsetI.length < maxAntecedentSize){
						rulesForRecursion.rules.add(newRule);
					}
				}
			}
			
			if(rulesForRecursion.rules.size() >1  ) {
				expandLeft(rulesForRecursion);
			}
		}
		// check the memory usage
    	MemoryLogger.getInstance().checkMemory();
	}


	private void expandRight(RightEquivalenceClass eclass, boolean firstTime) {
//		System.out.println("EXPANDRIGHT");
//		System.out.println("Trying to expand equivalence class " + eclass.toString());
		for(int w=0; w < eclass.rules.size()-1; w++){ // IMPORTANT : SIZE -1 BECAUSE THE LAST ONE HAS NOTHING LEFT FOR COMPARISON
			RightRule rule1 = eclass.rules.get(w);
			int d = rule1.itemsetJ[rule1.itemsetJ.length -1];
			RightEquivalenceClass rulesForRecursion= new RightEquivalenceClass(eclass.itemsetI, eclass.tidsI, eclass.occurencesI);

			// for each rule J != I
			for(int m=w+1; m < eclass.rules.size(); m++)	{
				RightRule rule2 = eclass.rules.get(m);
				
				int c = rule2.itemsetJ[rule2.itemsetJ.length -1];

				if(matrix.getCount(c, d) < minsupAbsolute) {
					candidatePrunedCount++;
					totalCandidateCount++;
					continue;
				}
				totalCandidateCount++;
				
				// CALCULATE TIDS OF  I ==> JC
				Set<Integer> tidsI_JC = new HashSet<Integer>();
				
	 			// for each sequence containing I
				Map<Integer, Occurence> mapC = mapItemCount.get(c);

				// depending of the relative size of tids of IUJ and tids of {c} we have two ways to calculate it
				if(rule1.tidsIJ.size() < mapC.size()) {
					// EARLY SKIP OPTIMIZATION:
					int remains = rule1.tidsIJ.size();
					for(Integer tid: rule1.tidsIJ){
						// Get the first and last occurences of C in that sequence
						Occurence occurenceC = mapC.get(tid);
						// if there is an occurence
			    		if(occurenceC != null) {
							Occurence occurenceI = eclass.occurencesI.get(tid);
			    			if(occurenceC.lastItemset > occurenceI.firstItemset){
			    				// add the tid of the sequence to the tidset of JU{c}
			    				tidsI_JC.add(tid);
			    			}
			    		}
			    		remains--;
			    		if(tidsI_JC.size() + remains < minsupAbsolute) {
			    			break;
			    		}
			    	}
				}else {
					// EARLY SKIP OPTIMIZATION:
					int remains = mapC.size();
					for(Entry<Integer, Occurence> entryC: mapC.entrySet()){
						int tid = entryC.getKey();
						// if there is an occurence
			    		if(rule1.tidsIJ.contains(tid)) {
							Occurence occurenceC = entryC.getValue();
							Occurence occurenceI = eclass.occurencesI.get(tid);
			    			if(occurenceC.lastItemset > occurenceI.firstItemset){
			    				// add the tid of the sequence to the tidset of JU{c}
			    				tidsI_JC.add(tid);
			    			}
			    		}
			    	}
					remains--;
		    		if(tidsI_JC.size() + remains < minsupAbsolute) {
		    			break;
		    		}
				}

				ArrayList<Double> minDistsIC_J = getMinDists(rule1.itemsetJ, eclass.itemsetI, tidsI_JC, c,"right");
				double decayedSupp = .0;
				for(Double dist : minDistsIC_J) decayedSupp += dist;
				
				// if the support of I ==> JU{c} is enough 
				//if(tidsI_JC.size() >= minsupAbsolute){
				if(decayedSupp >= minsupAbsolute){

					// CALCULATE THE OCCURENCES OF JU{c}
	    			Set<Integer> tidsJC = new HashSet<Integer>(rule1.tidsJ.size());
	    			Map<Integer, Occurence> occurencesJC = new HashMap<Integer, Occurence>();
	    			
	    			// for each sequence containing J
	    			if(rule1.tidsJ.size() < mapC.size()) {
		    			for(Integer tid: rule1.tidsJ){
		    				// Get the first and last occurences of C in that sequence
		    				Occurence occurrenceC = mapC.get(tid);
		    				// if there is an occurence
		    	    		if(occurrenceC != null){
		    	    			// add the tid of the sequence to the tidset of JU{c}
		    	    			tidsJC.add(tid);
		    	    			// calculate last occurence of JU{c} depending on if
		    	    			// the last occurence of J is before the last occurence
		    	    			// of c or not.
		    	    			Occurence occurenceJ = rule1.occurencesJ.get(tid);
		    	    			if(occurrenceC.lastItemset < occurenceJ.lastItemset){
		    	    				occurencesJC.put(tid, occurrenceC);
		    	    			}else{
		    	    				occurencesJC.put(tid, occurenceJ);
		    	    			}
		    	    		}
		    	    	}
	    			}else {
	    				for(Entry<Integer, Occurence> entryC: mapC.entrySet()){
							int tid = entryC.getKey();
							// if there is an occurence
				    		if(rule1.tidsJ.contains(tid)) {
								// add the tid of the sequence to the tidset of JU{c}
		    	    			tidsJC.add(tid);
		    	    			// calculate last occurence of JU{c} depending on if
		    	    			// the last occurence of J is before the last occurence
		    	    			// of c or not.
								Occurence occurrenceC = entryC.getValue();
		    	    			Occurence occurenceJ = rule1.occurencesJ.get(tid);
		    	    			if(occurrenceC.lastItemset < occurenceJ.lastItemset){
		    	    				occurencesJC.put(tid, occurrenceC);
		    	    			}else{
		    	    				occurencesJC.put(tid, occurenceJ);
		    	    			}
				    		}
	    				}
	    			}
	    			
	    			// Create rule I ==> J U{c} and calculate its confidence   
	    	    	// defined as:  sup(I -->J U{c}) /  sup(I)	
	    			double confI_JC = ((double)tidsI_JC.size()) / eclass.tidsI.size();
					int[] itemsetJC = new int[rule1.itemsetJ.length+1];
					System.arraycopy(rule1.itemsetJ, 0, itemsetJC, 0, rule1.itemsetJ.length);
					itemsetJC[rule1.itemsetJ.length]= c;
				
					// if the confidence is enough
					if(confI_JC >= minConfidence){
						// then it is a valid rule so save it
						saveRule(decayedSupp, tidsI_JC, confI_JC, eclass.itemsetI, itemsetJC);
					}
					// recursively try to expand the left and right side
					// of the rule
					RightRule rightRule =
							new RightRule(itemsetJC, tidsJC, tidsI_JC, occurencesJC);
					
					if(rightRule.itemsetJ.length < maxConsequentSize){
						rulesForRecursion.rules.add(rightRule);
					}
					
					if(eclass.itemsetI.length < maxAntecedentSize){
						LeftRule leftRule = new LeftRule(eclass.itemsetI, eclass.tidsI, tidsI_JC);
						store.register(leftRule, itemsetJC, tidsJC, eclass.occurencesI, occurencesJC); // register for left expansion
					}
				}
			}

			if(rulesForRecursion.rules.size() >1) {
				expandRight(rulesForRecursion, false);
			}
		}
    	// check the memory usage
    	MemoryLogger.getInstance().checkMemory();
    	
	}

	/**
	 * Get a list of minimum distances of each occurence of I U {C} ==> J i.e. I ==> J U {C} in a sequence
	 * @param tids the tids containing the rule I U {C} ==> J or I ==> J U {C} respectively
	 * @param itemsetI the left part of the rule
	 * @param itemsetJ the right part of the rule
	 * @param itemC the item, which extends the antecedent and consequent respectively
	 */
	private ArrayList<Double> getMinDists(int[] itemsetI, int[] itemsetJ, Set<Integer> tids, Integer itemC, String expand) {
		ArrayList<Double> minDists = new ArrayList<>();
		for (Integer tid : tids) {
			int minDist = Integer.MAX_VALUE;

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

			for (int item : itemsetJ) {
				occsJ.addAll(mapItemCount.get(item).get(tid).allItemsets);
			}
			Collections.sort(occsI);
			Collections.sort(occsC);
			Collections.sort(occsJ);
			//System.out.println("Occs I: " + occsI);
			//System.out.println("Occs C: " + occsC);
			//System.out.println("Occs J: " + occsJ);

			// if antecedent side gets expanded
			if (expand.equals("left")) {
				for (Short occJ : occsJ) {
					int distIC_J = Integer.MAX_VALUE;
					// Check minimum distance between I and J
					int minDistI = Integer.MAX_VALUE;
					for (Short x : occsI) {
						int temp = occJ - x;
						if (temp >= 0 && temp < minDistI) {
							minDistI = temp;
						}
					}
					// Check minimum distance between C and J
					int minDistC = Integer.MAX_VALUE;
					for (Short x : occsC) {
						int temp = occJ - x;
						if (temp >= 0 && temp < minDistC) {
							minDistC = temp;
						}
					}
					// If I and C both occur before J take the minimum distance of both
					if (minDistC < Integer.MAX_VALUE && minDistI < Integer.MAX_VALUE) {
						if (minDistC < minDistI) {
							distIC_J = minDistC;
						} else {
							distIC_J = minDistI;
						}
					}
					// Search for minimum distance between IUc and J over the whole sequence
					if (distIC_J < minDist) {
						minDist = distIC_J;
					}
				}
				// Since I and c both have to occur at least once before J we have the minimum gapsize
				// over all the sequences where IUc-->J occurs
				minDists.add(Math.pow(decayFactor, minDist));
			} else { // if consequent side gets expanded
				for (Short occI : occsI) {
					int distI_JC = Integer.MAX_VALUE;
					// Check minimum distance between I and J
					int minDistJ = Integer.MAX_VALUE;
					for (Short x : occsJ) {
						int temp = x - occI;
						if (temp >= 0 && temp < minDistJ) {
							minDistJ = temp;
						}
					}
					// Check minimum distance between C and J
					int minDistC = Integer.MAX_VALUE;
					for (Short x : occsC) {
						int temp = x - occI;
						if (temp >= 0 && temp < minDistC) {
							minDistC = temp;
						}
					}
					// If I and C both occur before J take the minimum distance of both
					if (minDistC < Integer.MAX_VALUE && minDistJ < Integer.MAX_VALUE) {
						if (minDistC < minDistJ) {
							distI_JC = minDistC;
						} else {
							distI_JC = minDistJ;
						}
					}
					// Search for minimum distance between I and JUc over the whole sequence
					if (distI_JC < minDist) {
						minDist = distI_JC;
					}
				}
				// Since I and c both have to occur at least once before J we have the minimum gapsize
				// over all the sequences where IUc-->J occurs
				minDists.add(Math.pow(decayFactor, minDist));
			}
		}
		return minDists;
	}


	/**
	 * This method calculate the frequency of each item in one database pass.
	 * Then it remove all items that are not frequent.
	 * @param database : a sequence database 
	 * @return A map such that key = item
	 *                         value = a map  where a key = tid  and a value = Occurence
	 * This map allows knowing the frequency of each item and their first and last occurence in each sequence.
	 */
	private Map<Integer, Map<Integer, Occurence>> calculateFrequencyOfEachItem(SequenceDatabase database) {
		// (1) Count the support of each item in the database in one database pass
		mapItemCount = new HashMap<Integer, Map<Integer, Occurence>>(); // <item, Map<tid, occurence>>
		
		// for each sequence in the database
		for(int k=0; k< database.size(); k++){
			Sequence sequence = database.getSequences().get(k);
			// for each itemset in that sequence
			for(short j=0; j< sequence.getItemsets().size(); j++){
				List<Integer> itemset = sequence.get(j);
				// for each item in that itemset
				for(Integer itemI : itemset){
					
					// get the map of occurences of that item
					Map<Integer, Occurence> occurences = mapItemCount.get(itemI);
					// if this map is null, create a new one
					if(occurences == null){
						occurences = new HashMap<Integer, Occurence>();
						mapItemCount.put(itemI, occurences);
						Occurence occurence = new Occurence(j, j);
						occurence.allItemsets.add(j);
						occurences.put(k, occurence);
					}else {
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
		}
		
		// return the map of occurences of items
		return mapItemCount;
	}


	private void generateMatrix(SequenceDatabase database) {
		// for each sequence
		for(Sequence sequence : database.getSequences()){
			
			// to remember which items have been processed
			Set<Integer> alreadyProcessed = new HashSet<Integer>();
			
			// for each itemset
			for(List<Integer> itemsetj : sequence.getItemsets()) {
				
				// for each item
				for(Integer itemk : itemsetj) {
					if(alreadyProcessed.contains(itemk) || mapItemCount.get(itemk).size() < minsupAbsolute){
						continue;
					}
					
					/// for item k we should update the matrix with each item co-occurring 
					Set<Integer> alreadyProcessedWithRespectToK = new HashSet<Integer>();
					for(List<Integer> itemsetjj : sequence.getItemsets()) {
						
						for(Integer     itemkk    : itemsetjj) {
							if(itemkk.equals(itemk) || alreadyProcessedWithRespectToK.contains(itemkk)
									||  mapItemCount.get(itemkk).size() < minsupAbsolute){
								continue;
							}
							
							matrix.increaseCountOfPair(itemk, itemkk);
							alreadyProcessedWithRespectToK.add(itemkk);
						}
					}
					// end second loop
					alreadyProcessed.add(itemk);
				}
			}
		}
	}
	


	/**
	 * Save a rule I ==> J to the output file
	 * @param tidsIJ the tids containing the rule
	 * @param confIJ the confidence
	 * @param itemsetI the left part of the rule
	 * @param itemsetJ the right part of the rule
	 * @throws IOException exception if error writing the file
	 */
	private void saveRule(double decayedSupIJ, Set<Integer> tidsIJ, double confIJ, int[] itemsetI, int[] itemsetJ) {
//		System.out.println("######### Saving new rule");
		// increase the number of rule found
		ruleCount++;
		
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

		try {
			writer.write(buffer.toString());
			writer.newLine();
//			System.out.println(buffer.toString());
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	public void setDecayFactor(double decayFactor) { this.decayFactor = decayFactor;}
	
	/**
	 * Print statistics about the last algorithm execution to System.out.
	 */
	public void printStats() {
		System.out.println("=============  ERMiner - STATS ========");
		System.out.println("Sequential rules count: " + ruleCount);
		System.out.println("Total time: " + (timeEnd - timeStart) + " ms");
		System.out.println("Candidates pruned (%)" + candidatePrunedCount + " of " + totalCandidateCount);
		System.out.println("Max memory: " + MemoryLogger.getInstance().getMaxMemory());
		System.out.println("==========================================");
	}


	/**
	 * Set the number of items that a rule antecedent should contain (optional).
	 * @param maxAntecedentSize the maximum number of items
	 */
	public void setMaxAntecedentSize(int maxAntecedentSize) {
		this.maxAntecedentSize = maxAntecedentSize;
	}


	/**
	 * Set the number of items that a rule consequent should contain (optional).
	 * @param maxConsequentSize the maximum number of items
	 */
	public void setMaxConsequentSize(int maxConsequentSize) {
		this.maxConsequentSize = maxConsequentSize;
	}

	public static void main (String[] args) {
		AlgoERMiner er = new AlgoERMiner();
		//read in all arguments:
		double sup = Double.parseDouble(args[0]);
		double conf = Double.parseDouble(args[1]);
		String input = args[2];
		String output = args[3];
		er.setMaxConsequentSize(Integer.parseInt(args[4]));
		er.setDecayFactor(Double.parseDouble(args[5]));
//		System.out.println("Received the following input");
//		for (String s  : args){
//			System.out.println(s);
//		}

		try {
			er.runAlgorithm(sup, conf, input, output);
			//er.runAlgorithm(input, output, sup, conf);
		} catch (Exception e) {
			System.out.println("ERROR in main: " + e);
		}
	}

}
