import java.util.*;

public class Original93 {
	static Random r;
	
	final byte seq[][];
	final int N, L, W;
	final int[] psp, psn;
	final int pspsum, psnsum;
	int[] best, cur;
	double best_f, cur_f;
	int sample_pos = 0;
	
	// positive model (correspond to cur)
	int[][] c;
	
	// negative model (correspond to cur)
	int[] d;
	int dsum;
	
	// constructor: initialization
	public Original93(Seqs s, int W, int[] pseudo_pos, int[] pseudo_neg) throws Exception {
		seq = s.seq;
		N = s.N; L = s.L;
		this.W = W;

		if (pseudo_pos == null) {
			pseudo_pos = new int[L];
		}
		if (pseudo_neg == null) {
			pseudo_neg = new int[L];
		}
		if (pseudo_pos.length != L) {
			throw new Exception("Invalid pseudocount");
		}
		if (pseudo_neg.length != L) {
			throw new Exception("Invalid pseudocount");
		}
		psp = pseudo_pos; psn = pseudo_neg;
		int psum = 0, nsum = 0;
		for (int i = 0; i < L; i++) {
			psum += psp[i];
			nsum += psn[i];
		}
		pspsum = psum; psnsum = nsum;
		
		for (int i = 0; i < N; i++) {
			if (seq[i].length < W) {
				throw new Exception("The input is infeasible: each sequence should not be shorter than W");
			}
		}
		
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
			double pj = (double) (d[j] + psn[j]) / (dsum + psnsum);
			for (int i = 0; i < W; i++) {
				double qij = (double) (c[i][j] + psp[j]) / (N + pspsum);
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
		
		sample_pos++; // next seq
		if (sample_pos == N) sample_pos = 0;
		cur_f = evaluate();
		if (cur_f > best_f) {
			for (int i = 0; i < N; i++) {
				best[i] = cur[i];
			}
			best_f = cur_f;
			System.out.println("BEST: " + best_f);
		}
	}
	
	// calculate A_x (fix N-1 variables and set the one variable to pos)
	// approximated (O(W))
	private double fix_and_evaluate(int pos) throws Exception {
		double sum = 0;
		for (int i = 0; i < W; i++) {
			int j = seq[sample_pos][pos + i];
			// Qustion: Should the pseudocount be discounted by (N-1)/N?
			double pj = (double) (d[j] + psn[j]) / (dsum + psnsum);
			double qij = (double) (c[i][j] + psp[j]) / (N - 1 + pspsum);
			if (pj == 0) {
				throw new Exception("zero occurence alphabet: must have pseudocount > 0");
			}
			sum += Math.log(qij / pj);
		}
		return Math.exp(sum);
	}
	
	public static void main(String[] args) throws Exception {
		r = new Random();

		int[] psp = new int[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		int[] psn = new int[] {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		Original93 gibbs = new Original93(Seqs.read("input"), 18, psp, psn);
		
		for (int i = 0; i < 40000; i++) {
			gibbs.sample();
		}
		
		for (int i = 0; i < gibbs.N; i++) {
			System.out.print(gibbs.best[i] + " ");
		}
		System.out.println();
	}
}
