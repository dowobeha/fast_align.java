// Copyright 2013 by Chris Dyer
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Ported to Java by Lane Schwartz
//
package edu.cmu.clab;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

public class Dict {

	public Dict() {
		b0_= "<bad0>";
		d_ = new HashMap<String,Integer>();
		words_ = new ArrayList<String>(1000);
	}

	public int max() { return words_.size(); }

	public static boolean is_ws(char x) {
		return (x == ' ' || x == '\t');
	}

	public void ConvertWhitespaceDelimitedLine(final String line, ArrayList<Integer> out) {
		int cur = 0;
		int last = 0;
		int state = 0;
		out.clear();
		while(cur < line.length()) {
			if (is_ws(line.charAt(cur++))) {
				if (state == 0) continue;
				out.add(Convert(line.substring(last, last + (cur - last - 1))));
				state = 0;
			} else {
				if (state == 1) continue;
				last = cur - 1;
				state = 1;
			}
		}
		if (state == 1)
			out.add(Convert(line.substring(last, last + (cur - last))));
	}

	public int Convert(final String word) {
		return Convert(word, false);
	}

	public int Convert(final String word, boolean frozen) {
		Integer i = d_.get(word);
		if (i == null) {
			if (frozen)
				return 0;
			words_.add(word);
			d_.put(word, words_.size());
			return words_.size();
		} else {
			return i;
		}
	}

	public String Convert(final int id) {
		if (id == 0) return b0_;
		return words_.get(id-1);
	}

	private String b0_;
	private ArrayList<String> words_;
	private Map<String,Integer> d_;

	public static void ReadFromFile(final String filename,
			Dict d,
			ArrayList<ArrayList<Integer>> src,
			Set<Integer> src_vocab) throws FileNotFoundException {
		src.clear();
		System.err.println("Reading from " + filename);
		Scanner in = new Scanner(new File(filename), "UTF-8");
		assert(in.hasNextLine());
		String line;
//		int lc = 0;
		while(in.hasNextLine()) {
			line = in.nextLine();
//			++lc;
			ArrayList<Integer> back = new ArrayList<Integer>();
			src.add(back);
			d.ConvertWhitespaceDelimitedLine(line, back);
			for (int i = 0; i < back.size(); ++i) src_vocab.add(back.get(i));
		}
	}

}
