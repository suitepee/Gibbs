import java.util.*;
import java.io.*;

public class Seqs {
	final int N; // sequence
	final int L; // alphabet
	byte seq[][];
	
	private Seqs(int N, int L) {
		this.N = N;
		this.L = L;
	}
	
	public static Seqs read(String file) throws Exception {
		Scanner sc = new Scanner(new File(file));
		int N = sc.nextInt(), L = sc.nextInt();
		Seqs s = new Seqs(N, L);
		s.seq = new byte[N][];
		HashMap<Character, Byte> residue = new HashMap<Character, Byte>();
		
		String tmp;
		while (sc.hasNext() && ((tmp = sc.nextLine()).length() == 0 || tmp.charAt(0) != '>'));
		for (int i = 0; i < N; i++) {
			String seq = "";
			while (sc.hasNext() && (tmp = sc.nextLine()).charAt(0) != '>') {
				seq += tmp;
			}
			s.seq[i] = new byte[seq.length()];
			for (int j = 0; j < seq.length(); j++) {
				Byte rid = residue.get(seq.charAt(j));
				if (rid == null) {
					rid = (byte) residue.size();
					residue.put(seq.charAt(j), rid);
				}
				s.seq[i][j] = rid;
			}
		}
		sc.close();
		return s;
	}
}
