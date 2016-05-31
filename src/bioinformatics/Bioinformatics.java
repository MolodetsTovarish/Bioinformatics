/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformatics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Misha
 */
public class Bioinformatics {

    // Returns the k-mers of k length found in the genome, and how often it appears
    public static Hashtable kMerCount(String genome, int k) {
        Hashtable<String, Integer> kMer_counts = new Hashtable();
        for (int i = 0; i < genome.length() - (k - 1); i++) {
            String k_mer = genome.substring(i, k + i);

            Integer currentCount = kMer_counts.get(k_mer);
            if (currentCount != null) {

                kMer_counts.put(k_mer, currentCount + 1);
            } else {
                kMer_counts.put(k_mer, 1);
            }
            System.out.println(k_mer);
        }
        System.out.println(kMer_counts.toString());
        return kMer_counts;
    }

    // Counts the number of times each kmer occurs in the hashtable, and adds that kmer's compliment to the total count
    public static Hashtable kMerUnifiedCount(String genome, int k) {
        Hashtable<String, Integer> kMer_unified_counts = new Hashtable();
        String currentKey;
        Integer currentCount = 0;
        Hashtable<String, Integer> kMer_Counts = kMerCount(genome, k);
        Set<String> keys = kMer_Counts.keySet();

        //Obtaining iterator over set entries
        Iterator<String> itr = keys.iterator();

        //Displaying Key and value pairs
        while (itr.hasNext()) {
            // Getting Key
            currentKey = itr.next();
            currentCount = kMer_Counts.get(currentKey);

            if (!kMer_unified_counts.containsKey(currentKey)) {
                String transcriptedCurrentKey = transcriptGenome(currentKey);
                Integer transcriptedCurrentValue = kMer_unified_counts.get(transcriptedCurrentKey);
                if (transcriptedCurrentValue != null) {
                    kMer_unified_counts.put(transcriptedCurrentKey, transcriptedCurrentValue + currentCount);
                } else {
                    kMer_unified_counts.put(currentKey, currentCount);
                }
            }
            //System.out.println("Key: " + currentKey + ", Value: " + kMer_unified_counts.get(currentKey));
        }

        return kMer_unified_counts;

    }

    // Takes the top N (N being how many results you want to see) k-mers from the genome, from most to least frequent
    public static String topN(Hashtable kMer_count, int n) {
        ArrayList kMer_list = sortValue(kMer_count);
        return kMer_list.subList(0, n).toString();
    }

    public static ArrayList sortValue(Hashtable<?, Integer> t) {

        //Transfer as List and sort it
        ArrayList<Map.Entry<?, Integer>> l = new ArrayList(t.entrySet());
        Collections.sort(l, new Comparator<Map.Entry<?, Integer>>() {

            public int compare(Map.Entry<?, Integer> o1, Map.Entry<?, Integer> o2) {
                return o2.getValue().compareTo(o1.getValue());
            }
        });

        System.out.println(l);
        return l;
    }

    public static String topKey(Hashtable<String, Integer> kMer_count) {
        Set<String> keys = kMer_count.keySet();
        //Obtaining iterator over set entries
        Iterator<String> itr = keys.iterator();
        Integer maxValue = 0;
        Integer currentValue = 0;
        String topKey = null;
        String currentKey = null;

        //Displaying Key and value pairs
        while (itr.hasNext()) {
            // Getting Key
            currentKey = itr.next();
            currentValue = kMer_count.get(currentKey);

            if (currentValue > maxValue) {
                maxValue = currentValue;
                topKey = currentKey + " - " + maxValue.toString();
            }

            System.out.println("Key: " + currentKey + ", Value: " + kMer_count.get(currentKey));
        }

        return topKey;
    }

    public static String transcriptGenome(String genome) {
        String complimentaryString = "";
        for (int i = genome.length() - 1; i >= 0; i--) {
            complimentaryString = complimentaryString + transcriptionRules(genome.charAt(i));
        }
        return complimentaryString;
    }

    public static char transcriptionRules(char nucleotide) {
        if (nucleotide == 'a') {
            return 't';
        } else if (nucleotide == 't') {
            return 'a';
        } else if (nucleotide == 'c') {
            return 'g';
        } else if (nucleotide == 'g') {
            return 'c';
        } else {
            //return Character.MIN_VALUE;
            throw new IllegalArgumentException("Invalid symbol");
        }
    }
    
    public static ArrayList findKMerPositions(String genome, String kMer){
        ArrayList<Integer> kMer_positions = new ArrayList<Integer>();
        int kMer_length = kMer.length();
        for (int i = 0; i < genome.length() - (kMer_length - 1); i++) {
            if (genome.substring(i, i + kMer_length).equals(kMer)){ 
                kMer_positions.add(i);
            }
        }
        return kMer_positions;
    }
    
    /*public static ArrayList findKMerClumps(String genome, int k, int l) {
        //ArrayList<Integer> kMer_positions = findKMerPositions(genome, kMer);
        //ArrayList kMer_list = new ArrayList();
        //ArrayList clump_positions = new ArrayList();
        Hashtable<String, Integer> clump_count = new Hashtable();
        String currentKey;
        Integer currentCount = 0;
        int genome_size = genome.length();
        
        for (int i = 0; i < genome.length() - (k - 1); i++) {
            String current_kMer = genome.substring(i, k + i);
            currentCount = clump_count.get(current_kMer);
            
            if (!kMer_list.contains(current_kMer)){
                kMer_list.add(current_kMer);
                
                if ()
            }

            
        }
    }
    */
    
    public static ArrayList findKMerClumps(String genome, int k, int l){
        Hashtable<String, Integer> kMer_positions = new Hashtable();
        ArrayList clump_list = new ArrayList();
        int genome_size = genome.length();
        int previous_positions;
        
        //loops until it reaches te end of the genome minus the kmer length
        for (int i = 0; i < genome_size - (k - 1); i++) {
            String current_kMer = genome.substring(i, k + i);
            
            // if the clump list has the current kmer...
            if (!clump_list.contains(current_kMer)){
                // if the kmer positions list doesn't have the current kmer...
                if (!kMer_positions.contains(current_kMer)){
                    // put that kmer and its position in the kmer positions list
                    kMer_positions.put(current_kMer, i);
                }
                // otherwise, if the kmer list doesn't have current kmer...
                else{
                
                    // if the current index minus the index of the previous position of kmer is less than L...
                    if ((i - kMer_positions.get(current_kMer)) <= l){
                        // add the kmer to the clump list
                        clump_list.add(current_kMer);
                        
                    }
                    // otherwise, replace the previous position of kmer with the current one
                    else{
                        kMer_positions.replace(current_kMer, i);
                    }
                //}
                }
            }
        }
        
        return clump_list;
    }
    
    public static void main(String[] args) {
       // System.out.print(topKey(kMerCount("acaaccccac", 2)));
        //System.out.println(kMerUnifiedCount("ccgg", 2));
        //System.out.println(topKey(kMerUnifiedCount("ccgg", 2)));
        //System.out.println("\n" + transcriptGenome("att"));
        
        //System.out.println(findKMerPositions("acttgactattt", "tt"));
        
        System.out.println(findKMerClumps("aataaaaaat", 2, 2));
    }

}
