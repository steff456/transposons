package ngsep.genome;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.QualifiedSequence;

public class Transposon {
	
	private Distribution distrHits = new Distribution(0, 100, 1);
	
	private ReferenceGenome genome;
	
	private ReferenceGenomeFMIndex fm;
	
	private DefaultKmersMapImpl kmersMap;
	
	private HashMap<String, List<GenomicRegionImpl>> candidates;
		
	public void run() {
		// Kmers - Subsequence
		int numSequences = genome.getNumSequences();
		for (int i = 0; i < numSequences; i++) {
			QualifiedSequence qs = genome.getSequenceByIndex(i);
			CharSequence seq = qs.getCharacters();
			processSequence(seq, qs.getName(), fm);			
		}
		//TODO: process tandem repeats
		//TODO: Find candidates that are nearer than 5kb
		//TODO: Find LTR form
	}
	
	public void printResults() {
		distrHits.printDistribution(System.out);		
	}
	
	public void processSequence(CharSequence seq, String name, ReferenceGenomeFMIndex fm) {
		System.out.printf("Processing Sequence %s \n", name);
		//Subsequence 20bp
		int lengthKmer = 20;
		for (int i = 0; i + lengthKmer < seq.length(); i+=10) {
			CharSequence kmer = seq.subSequence(i, (i+lengthKmer));
			List<ReadAlignment> hits = fm.search(kmer.toString());
			distrHits.processDatapoint(hits.size());
			if(hits.size() > 100 ) {
				// If the kmer is more than x times
				kmersMap.addOcurrance(kmer);
				GenomicRegionImpl temp = new GenomicRegionImpl(kmer.toString(), i, i+lengthKmer);
				if(!candidates.keySet().contains(kmer.toString())) {
					List n = new ArrayList<GenomicRegionImpl>();
					n.add(temp);
					candidates.put(kmer.toString(), n);
				}
				else {
					List n = candidates.get(kmer.toString());
					n.add(temp);
					candidates.put(kmer.toString(), n);
				}
			}
		}
	}
	
	public static void main(String[] args) throws Exception {
		Transposon instance = new Transposon();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// Load short tandem repeats from the .fa file
		//TODO: Write the loader for short tandem repeats
		
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		
		instance.kmersMap = new DefaultKmersMapImpl();
		
		instance.candidates = new HashMap();
		
		instance.run();
		instance.printResults();
		System.out.println(instance.kmersMap.size());
	}
}
