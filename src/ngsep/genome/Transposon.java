package ngsep.genome;
import java.io.IOException;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.QualifiedSequence;

public class Transposon {
	
	private Distribution distrHits = new Distribution(0, 100, 1);
	
	private ReferenceGenome genome;
	
	private ReferenceGenomeFMIndex fm;
	
	public void run() {
		// Kmers - Subsequence
		int numSequences = genome.getNumSequences();
		for (int i = 0; i < numSequences; i++) {
			QualifiedSequence qs = genome.getSequenceByIndex(i);
			CharSequence seq = qs.getCharacters();
			processSequence(seq, qs.getName(), fm);
		}
		
	}
	
	public void printResults() {
		distrHits.printDistribution(System.out);		
	}
	
	public void processSequence(CharSequence seq, String name, ReferenceGenomeFMIndex fm) {
		//Subsequence 20bp
		System.out.printf("Processing Sequence %s \n", name);
		int lengthKmer = 20;
		for (int i = 0; i + lengthKmer < seq.length(); i+=10) {
			CharSequence kmer = seq.subSequence(i, (i+lengthKmer));
			List<ReadAlignment> hits = fm.search(kmer.toString());
			distrHits.processDatapoint(hits.size());
			if(hits.size() > 100 )
				System.out.println("Seq: " + kmer + " Hits: " + hits.size());
		}
	}
	
	public static void main(String[] args) throws Exception {
		Transposon instance = new Transposon();
		// Load the genome from a .fa file
		instance.genome = new ReferenceGenome(args[0]);
		// FM Index
		instance.fm = new ReferenceGenomeFMIndex(instance.genome);
		
		instance.run();
		instance.printResults();
	}
}
