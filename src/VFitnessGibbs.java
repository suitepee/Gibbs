import java.util.*;

public class VFitnessGibbs {
	static Random r;
	
	// variables from input
	final byte seq[][];
	final int N, L, W;
	
	// the weight of pseudocount
	final int psp, psn;
	
	// solution info
	int[] best, cur;
	double best_f, cur_f;

	// positive model from counting (correspond to cur)
	int[][] c; // c[W][L]
	
	// negative model from counting (correspond to cur)
	int[] d; // d[L]
	int dsum; // d[0] + ... + d[L-1]
	
	// gibbs sampler position [0...N)
	int sample_pos = 0;
	
	// constants for VFitpolicy (tuning needed)
	final double crmax, crmin;
	final double crunit;
	
	// variables for VFitpolicy
	double count_ratio;
	int loc = 0, all = 0;
	boolean updated = false;
	
	// constructor (initialize)
	public VFitnessGibbs(Seqs s, int W, int pseudo_pos, int pseudo_neg) throws Exception {
		seq = s.seq;
		N = s.N; L = s.L;
		this.W = W;

		psp = pseudo_pos; psn = pseudo_neg;
		
		for (int i = 0; i < N; i++) {
			if (seq[i].length < W) {
				throw new Exception("The input is infeasible: each sequence should not be shorter than W");
			}
		}
		
		crmax = N * N;
		crmin = (double) L / N;
		crunit = (crmax - crmin) / N;
		count_ratio = crmax;
		
		best = new int[N];
		cur = new int[N];
		for (int i = 0; i < N; i++) {
			int rint = r.nextInt(seq[i].length - W + 1);
			best[i] = rint;
			cur[i] = rint;
		}
		
		c = new int[W][L];
		d = new int[L];
		dsum = 0;
		for (int i = 0; i < N; i++) {
			// count positive
			byte[] cs = seq[i];
			for (int j = 0; j < W; j++) {
				c[j][cs[cur[i] + j]]++;
			}
			
			// count negative
			for (int j = 0; j < cur[i]; j++) {
				d[cs[j]]++;
				dsum++;
			}
			for (int j = cur[i] + W; j < seq[i].length; j++) {
				d[cs[j]]++;
				dsum++;
			}
		}
		best_f = evaluate();
		cur_f = best_f;
		System.out.println("INIT: " + best_f);
	}
	
	// calculate cur's fitness in O(W*L)
	private double evaluate() throws Exception {
		double sum = 0;
		for (int j = 0; j < L; j++) {
			double pj = (double) (d[j] + psn) / (dsum + psn * L);
			for (int i = 0; i < W; i++) {
				double qij = (double) (c[i][j] + psp) / (N + psp * L);
				if (pj == 0) {
					throw new Exception("zero occurence alphabet: must have pseudocount > 0");
				}
				if (c[i][j] != 0) {
					sum += c[i][j] * Math.log(qij / pj);
				}
			}
		}
		return sum;
	}
	
	public void sample() throws Exception {
		byte[] cs = seq[sample_pos];
		// subtract current position from the models
		for (int j = 0; j < W; j++) {
			c[j][cs[cur[sample_pos] + j]]--;
		}
		for (int j = 0; j < cur[sample_pos]; j++) {
			d[cs[j]]--;
			dsum--;
		}
		for (int j = cur[sample_pos] + W; j < cs.length; j++) {
			d[cs[j]]--;
			dsum--;
		}
		
		// sampling
		int range = cs.length - W + 1; // pos is in [0..range)
		double[] ratio = new double[range];
		double rsum = 0;
		// O(range*W)
		for (int i = 0; i < range; i++) {
			ratio[i] = fix_and_evaluate(i);
			rsum += ratio[i];
		}
		int smlnum = 0;
		for (int i = 0; i < range; i++) {
			if (ratio[i] <= ratio[cur[sample_pos]]) smlnum++;
		}
		
		double rand = r.nextDouble() * rsum;
		int sel;
		for (sel = 0; sel < range - 1; sel++) {
			if (rand < ratio[sel]) break;
			else {
				rand -= ratio[sel];
			}
		}
		
		cur[sample_pos] = sel;
		// apply new position to the models
		for (int j = 0; j < W; j++) {
			c[j][cs[cur[sample_pos] + j]]++;
		}
		for (int j = 0; j < cur[sample_pos]; j++) {
			d[cs[j]]++;
			dsum++;
		}
		for (int j = cur[sample_pos] + W; j < cs.length; j++) {
			d[cs[j]]++;
			dsum++;
		}
		
		if (range == smlnum) {
			loc++;
		}
		all++;
		
		// evaluate fitness and update if needed
		cur_f = evaluate();
		if (cur_f > best_f) {
			for (int i = 0; i < N; i++) {
				best[i] = cur[i];
			}
			best_f = cur_f;
			updated = true;
			System.out.println("BEST: " + best_f);
		}

		// policy
		sample_pos++; // next seq
		if (sample_pos == N) {
			sample_pos = 0;
			if (!updated) {
				if (all == loc) {
					count_ratio = crmax;
				} else {
					count_ratio -= crunit;
					if (count_ratio < crmin) count_ratio = crmin;
				}
			}
			updated = false;
			all = 0;
			loc = 0;
		}
	}
	
	// calculate A_x (fix N-1 variables and set the one variable to pos)
	// approximated (O(W))
	private double fix_and_evaluate(int pos) throws Exception {
		double sum = 0;
		for (int i = 0; i < W; i++) {
			int j = seq[sample_pos][pos + i];
			// Question: Should the pseudocount be discounted by (N-1)/N?
			double pj = (double) (d[j] + psn) / (dsum + psn * L);
			double qij = (double) (c[i][j] + count_ratio * N / L) / (N + count_ratio * N);
			if (pj == 0) {
				throw new Exception("zero occurence alphabet: must have pseudocount > 0");
			}
			sum += Math.log(qij / pj);
		}
		return Math.exp(sum);
	}
	
	public static void main(String[] args) throws Exception {
		long seed = System.currentTimeMillis();
		r = new Random(seed);
		VFitnessGibbs gibbs = new VFitnessGibbs(Seqs.read("input"), 18, 1, 1);
		
		long t = System.currentTimeMillis();
		while (true) {
			for (int i = 0; i < 3000; i++) {
				gibbs.sample();
			}
			if (System.currentTimeMillis() - t > 30000) break;
		}
		
		for (int i = 0; i < gibbs.N; i++) {
			System.out.print(gibbs.best[i] + " ");
		}
		System.out.println();
		System.out.println("random seed: " + seed);
	}
}
