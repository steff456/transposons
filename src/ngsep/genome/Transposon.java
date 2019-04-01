package ngsep.genome;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.QualifiedSequence;

public class Transposon {
	
	public class Tuple<X, Y> { 
		  public final X x; 
		  public final Y y; 
		  public Tuple(X x, Y y) { 
		    this.x = x; 
		    this.y = y; 
		  } 
	}
		
	private Distribution distrHits = new Distribution(0, 100, 1);
	
	private ReferenceGenome genome;
	
	private ReferenceGenomeFMIndex fm;
	
	private DefaultKmersMapImpl kmersMap;
	
	private HashMap<String, List<GenomicRegionImpl>> candidates;
	
	private List<Tuple<GenomicRegionImpl, GenomicRegionImpl>> ltr;
	
	private int count = 0;
		
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
		findNearCandidates();
		//TODO: Find LTR form
	}
	
	public void findNearCandidates() {
		for(String candidate: candidates.keySet()) {
			List<GenomicRegionImpl> genomicRegions = candidates.get(candidate);
			// Compare the distance among the genomic regions
			for (int i = 0; i < genomicRegions.size(); i++) {
				GenomicRegionImpl actRegion = genomicRegions.get(i);
				for (int j = i+1; j < genomicRegions.size(); j++) {
					GenomicRegionImpl secRegion = genomicRegions.get(j);
					int distance = GenomicRegionPositionComparator.getInstance().compare(actRegion, secRegion);
					if (distance >= 5000) {
						ltr.add(new Tuple<GenomicRegionImpl, GenomicRegionImpl>(actRegion, secRegion));
					}
				}
			}
			
		}
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
			if(hits.size() > 10 ) {
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
				count++;
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
		
		instance.ltr = new ArrayList<>();
		
		instance.run();
		instance.printResults();
		System.out.println(instance.kmersMap.size());
		System.out.println(instance.count);
		System.out.println(instance.ltr.size());
	}
}
