package ngsep.genome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.math.Distribution;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.QualifiedSequence;

public class TransposonFinder {
	
	/**
	 * Genome of interest
	 */
	private ReferenceGenome genome;
	
	/**
	 * FMIndex of the genome of interest
	 */
	private ReferenceGenomeFMIndex fm;
		
	private Map<String, List<GenomicRegion>> STRs;
	
	private Distribution distrHits = new Distribution(0, 100, 1);
	
	private HashMap<String, List<GenomicRegionImpl>> candidates;
	
	public void run() {
		// Kmers per subsequence
		int numSequences = genome.getNumSequences();
		for (int i = 0; i < numSequences; i++) {
			QualifiedSequence qs = genome.getSequenceByIndex(i);
			CharSequence seq = qs.getCharacters();
			processSequence(seq, qs.getName(), fm);			
		}
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
			}
		}
	}

	public static void main(String[] args) throws IOException {
		TransposonFinder instance = new TransposonFinder();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// Load short tandem repeats from the .fa file
		SimpleGenomicRegionFileHandler loader = new SimpleGenomicRegionFileHandler();
		instance.STRs = loader.loadRegionsAsMap(args[1]);
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		
		instance.run();
	}
}