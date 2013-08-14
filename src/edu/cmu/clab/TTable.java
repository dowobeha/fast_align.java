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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;

public class TTable {

	private static class Md {
		static double digamma(double x) {
			double result = 0, xx, xx2, xx4;
			for ( ; x < 7; ++x)
				result -= 1/x;
			x -= 1.0/2.0;
			xx = 1.0/x;
			xx2 = xx*xx;
			xx4 = xx2*xx2;
			result += Math.log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
			return result;
		}		
	}
	public TTable() {}
	double prob(final int e, final int f) {
		Map<Integer,Double> cpd = ttable.get(e);
		if (cpd != null) {
			Double value = cpd.get(f);
			if (value == null) return 1e-9;
			return value;
		} else {
			return 1e-9;
		}
	}
	public void Increment(final int e, final int f) {
		Increment(e,f,1);
	}
	public void Increment(final int e, final int f, double x) {
		Map<Integer,Double> map = counts.get(e);
		if (map==null) {
			map = new HashMap<Integer,Double>();
			counts.put(e, map);
		}
		Double value = map.get(f);
		if (value==null) value = 0.0;
		value += x;
		map.put(f, value);
	}
	public void NormalizeVB(final double alpha) {
		Map<Integer,Map<Integer,Double>> tmp = ttable;
		ttable = counts;
		counts = tmp;
		for (Iterator<Map.Entry<Integer,Map<Integer,Double>>> ttableIt = ttable.entrySet().iterator();
				ttableIt.hasNext(); ) {
			Map.Entry<Integer,Map<Integer,Double>> cit = ttableIt.next();
			double tot = 0;
			Map<Integer,Double> cpd = cit.getValue();
			for (Iterator<Map.Entry<Integer, Double>> cpdIt = cpd.entrySet().iterator(); cpdIt.hasNext(); ) {
				Map.Entry<Integer, Double> it = cpdIt.next();
				tot += it.getValue() + alpha;
			}
			for (Iterator<Map.Entry<Integer, Double>> cpdIt = cpd.entrySet().iterator(); cpdIt.hasNext(); ) {
				Map.Entry<Integer, Double> it = cpdIt.next();	
				it.setValue(Math.exp(Md.digamma(it.getValue() + alpha) - Md.digamma(tot)));
			}

		}
		counts.clear();
	}
	void Normalize() {
		Map<Integer,Map<Integer,Double>> tmp = ttable;
		ttable = counts;
		counts = tmp;
		for (Iterator<Map.Entry<Integer,Map<Integer,Double>>> ttableIt = ttable.entrySet().iterator();
				ttableIt.hasNext(); ) {
			Map.Entry<Integer,Map<Integer,Double>> cit = ttableIt.next();
			double tot = 0;
			Map<Integer,Double> cpd = cit.getValue();
			for (Iterator<Map.Entry<Integer, Double>> cpdIt = cpd.entrySet().iterator(); cpdIt.hasNext(); ) {
				Map.Entry<Integer, Double> it = cpdIt.next();
				tot += it.getValue();
			}
			for (Iterator<Map.Entry<Integer, Double>> cpdIt = cpd.entrySet().iterator(); cpdIt.hasNext(); ) {
				Map.Entry<Integer, Double> it = cpdIt.next();	
				it.setValue(it.getValue() / tot);
			}			
		}
		counts.clear();
	}
	/** adds counts from another TTable - probabilities remain unchanged */
	public TTable add(final TTable rhs) {
		for (Iterator<Map.Entry<Integer,Map<Integer,Double>>> rhsCountsIt = rhs.counts.entrySet().iterator();
				rhsCountsIt.hasNext(); ) {
			Map.Entry<Integer,Map<Integer,Double>> it = rhsCountsIt.next();
			final Map<Integer,Double> cpd = it.getValue();
			Map<Integer,Double> tgt = counts.get(it.getKey());
			for (Iterator<Map.Entry<Integer, Double>> jIt = cpd.entrySet().iterator(); jIt.hasNext(); ) {
				Map.Entry<Integer, Double> j = jIt.next();
				tgt.put(j.getKey(), tgt.get(j.getKey()) + j.getValue());
			}
		}
		return this;
	}
	public void ExportToFile(final String filename, Dict d) {
		PrintStream file = null;
		try {
			file = new PrintStream(new FileOutputStream(filename), true, "UTF-8");
			for (Iterator<Map.Entry<Integer,Map<Integer,Double>>> ttableIt = ttable.entrySet().iterator();
					ttableIt.hasNext(); ) {
				Map.Entry<Integer,Map<Integer,Double>> cit = ttableIt.next();
				final String a = d.Convert(cit.getKey());
				Map<Integer,Double> cpd = cit.getValue();
				for (Iterator<Map.Entry<Integer, Double>> cpdIt = cpd.entrySet().iterator(); cpdIt.hasNext(); ) {
					Map.Entry<Integer, Double> it = cpdIt.next();
					final String b = d.Convert(it.getKey());
					Double c = it.getValue();
					file.println(a + '\t' + b + '\t' + c);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				file.close();
			} catch (Exception ex) {}
		}
	}
	boolean ImportFromFile(String filename, char delim, Dict d) {
		Scanner in = null;
		try {
			in = new Scanner(new File(filename));

			in.useDelimiter("[\t\n\r]+");
			if (! in.hasNext()) {
				return false;
			} else {
				while(in.hasNext()) {

					String sourceWord=null;
					String targetWord=null;
					String valueString = null;

					boolean success = true;
					success &= in.hasNext(); if (success) sourceWord=in.next();
					success &= in.hasNext(); if (success) targetWord=in.next();
					success &= in.hasNext(); if (success) valueString=in.next();

					if (success) {
						int source = d.Convert(sourceWord);
						int target = d.Convert(targetWord);
						double value = Double.valueOf(valueString);
						Map<Integer,Double> map = ttable.get(source);
						if (map==null) {
							map = new HashMap<Integer,Double>();
							ttable.put(source, map);
						}
						map.put(target, value);
					} else {
						return false;
					}
				}
			}
			return true;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return false;
		} finally {
			try {
				if (in!=null) { in.close(); }
			} catch (Exception e) {}
		}
	}
	public Map<Integer,Map<Integer,Double>> ttable = new HashMap<Integer,Map<Integer,Double>>();
	public Map<Integer,Map<Integer,Double>> counts = new HashMap<Integer,Map<Integer,Double>>();	
}
