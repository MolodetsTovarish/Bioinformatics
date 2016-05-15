/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bioinformatics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
    
    public static Hashtable<String, Integer> kMer_counts = new Hashtable();
    
    // Returns the k-mers of k length found in the genome, and how often it appears
    public static Hashtable kMerCount(String genome, int k) {
        for (int i=0; i < genome.length() - (k-1); i++){
            String k_mer = genome.substring(i, k+i);
            
            Integer currentCount = kMer_counts.get(k_mer);
            if (currentCount != null) {
         
                kMer_counts.put(k_mer, currentCount+1);
            }
            else
                kMer_counts.put(k_mer, 1);
                System.out.println(k_mer);
        }
        System.out.println(kMer_counts.toString());
        return kMer_counts;
    }
    
    // Takes the top N (N being how many results you want to see) k-mers from the genome, from most to least frequent
    public static String topN(Hashtable kMer_count, int n){
        ArrayList kMer_list = sortValue(kMer_count);
        return kMer_list.subList(0, n).toString();
    }
    
    public static ArrayList sortValue(Hashtable<?, Integer> t){

       //Transfer as List and sort it
       ArrayList<Map.Entry<?, Integer>> l = new ArrayList(t.entrySet());
       Collections.sort(l, new Comparator<Map.Entry<?, Integer>>(){
       
         public int compare(Map.Entry<?, Integer> o1, Map.Entry<?, Integer> o2) {
             return o2.getValue().compareTo(o1.getValue());
        }});
       
       System.out.println(l);
       return l;
    }
    
    public static String topKey(Hashtable<String, Integer> kMer_count){
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
           
           if (currentValue > maxValue){
               maxValue = currentValue;
               topKey = currentKey + " - " + maxValue.toString();
           }
           
           System.out.println("Key: "+ currentKey+", Value: "+kMer_count.get(currentKey));
        } 
    
        return topKey;
    }
    
    public static String transcriptGenome(String genome){
        String complimentaryString = "";
        for (int i = genome.length()-1; i >= 0; i--){
            complimentaryString = complimentaryString + transcriptionRules(genome.charAt(i));
        }
        return complimentaryString;
    }
    
    public static char transcriptionRules(char nucleotide){
        if (nucleotide == 'a'){
            return 't';
        }
        else if (nucleotide == 't'){
            return 'a';
        }
        
        else if (nucleotide == 'c'){
            return 'g';
        }
        
        else if (nucleotide == 'g'){
            return 'c';
        }
        else {
        //return Character.MIN_VALUE;
        throw new IllegalArgumentException("Invalid symbol");
        }
    }
            
    public static void main(String[] args) {
        System.out.print(topKey(kMerCount("acaaccccac", 2)));
        System.out.println("\n" + transcriptGenome("aaaccaatttgg"));
    }
    
}
