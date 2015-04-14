/*
 * testWM.cpp
 *
 *  Created on: 2015年4月5日
 *      Author: John
 */
#include<WaveletMatrix.h>
#include<Array.h>
#include <boost/algorithm/string.hpp>
using namespace cds_static;
const int NUM = 2;

size_t *mystarts = new size_t[NUM] { 0, 4 /*, 223*/};
size_t *myends = new size_t[NUM] { 3, 12 /*, 227*/};

int main(int argc, char** argv) {
	cout << "hello world" << endl;
	const char* invlist = argv[1];
	const char* invlistfreq = argv[2];
	const char* vocab = argv[3];
	const char* doclens_file = argv[4];

	ifstream wordsfile;
	wordsfile.open(vocab);
	vector<string> words;
	string line;

	//读取词典文件
	while (wordsfile.good()) {
		getline(wordsfile, line);
		words.push_back(line);
	}
	uint count = words.size();
	for (int i = 0; i < count; i++)
		cout << i << ":" << words[i] << endl;
	//最后一行为空白，需要弹出
	words.pop_back();
	wordsfile.close();

	//读取倒排链
	ifstream docfile;
	docfile.open(invlist);
	vector<vector<int> > result;
	vector<string> strs;
	int length = 0;
//	int times=0;
	while (docfile.good()) {
		vector<int> r;
		getline(docfile, line);
		boost::split(strs, line, boost::is_any_of(" "));
		int doc_size = atoi(strs[0].c_str());
		length += doc_size;
		for (int i = 1; i <= doc_size; i++)
			r.push_back(atoi(strs[i].c_str()));
//				cout << times++ << "r.size():" << r.size() << endl;
//				for (int i = 0; i < r.size(); i++)
//					cout << r[i] << endl;

		result.push_back(r);
	}
	//最后一行为空白，弹出最后一行
	result.pop_back();
	docfile.close();

	//生成BitSequence
	BitString *bs=new BitString(length);
	int m2a=0;
	for (int i = 0; i < words.size(); i++) {
		bs->setBit(m2a);
		m2a = m2a + result[i].size();
	}

	BitSequenceRG *bsrg = new BitSequenceRG(*bs, 2);
	ofstream file("test.txt");
	bsrg->save(file);


//	//生成WM并测试
//	uint *sequence = new uint[length];
//	size_t m2a = 0;
//	for (int i = 0; i < result.size(); i++) {
//		for (int j = 0; j < result[i].size(); j++) {
//			sequence[m2a] = result[i][j];
//			m2a++;
//		}
//	}
//
//	Array *A = new Array(sequence, length);
//	MapperNone * map = new MapperNone();
//	BitSequenceBuilder * bsb = new BitSequenceBuilderRG(2);
//	WaveletMatrix* seq = new WaveletMatrix(*A, bsb, map);
//
//	seq->n_range_intersect_aux(mystarts, myends, NUM);

}
