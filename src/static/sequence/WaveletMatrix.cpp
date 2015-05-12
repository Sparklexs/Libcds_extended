/* WaveletMatrix.cpp
 * Copyright (C) 2012, Francisco Claude & Gonzalo Navarro, all rights reserved.
 *
 * WaveletMatrix definition
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <WaveletMatrix.h>

namespace cds_static {

WaveletMatrix::WaveletMatrix(const Array &symbols2, BitSequenceBuilder * bmb,
		Mapper * am) :
		Sequence(n) {
	bmb->use();
	n = symbols2.getLength();
	uint *symbols = new uint[n];
	this->am = am;
	am->use();
	for (size_t i = 0; i < n; i++)
		symbols[i] = am->map(symbols2.getField(i));
	max_v = max_value(symbols, n);
	height = bits(max_v);

	OCC = new uint[max_v + 2];
	for (uint i = 0; i <= max_v + 1; i++)
		OCC[i] = 0;
	for (size_t i = 0; i < n; i++)
		OCC[symbols[i] + 1]++;

	size_t to_add = 0;
	for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0)
			to_add++;

	uint * new_symb = new uint[n + to_add];
	for (size_t i = 0; i < n; i++)
		new_symb[i] = symbols[i];
	delete[] symbols;

	to_add = 0;
	for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0) {
			OCC[i]++;
			new_symb[n + to_add] = i - 1;
			to_add++;
		}

	size_t new_n = n + to_add;
	for (uint i = 1; i <= max_v + 1; i++)
		OCC[i] += OCC[i - 1];
	this->n = new_n;

	uint **_bm = new uint*[height];
	for (uint i = 0; i < height; i++) {
		_bm[i] = new uint[new_n / W + 1];
		for (size_t j = 0; j < new_n / W + 1; j++)
			_bm[i][j] = 0;
	}

	build_level(_bm, new_symb, new_n, NULL);
	bitstring = new BitSequence*[height];
	//bitstring[0]是最高位，依次递减
	//bitstring存储了每层的bitmap，就可以知道该symbol在对应的bit位是0 or 1
	C = new uint[height];
	for (uint i = 0; i < height; i++) {
		bitstring[i] = bmb->build(_bm[i], new_n);
		C[i] = bitstring[i]->rank0(new_n - 1);
		delete[] _bm[i];
	}
	delete[] _bm;
	// delete [] oc;
	bmb->unuse();

	this->length = n;
}

WaveletMatrix::WaveletMatrix(uint * symbols, size_t n, BitSequenceBuilder * bmb,
		Mapper * am, bool deleteSymbols) :
		Sequence(n) {
	bmb->use();
	this->n = n;
	this->am = am;
	am->use();
	for (size_t i = 0; i < n; i++)
		symbols[i] = am->map(symbols[i]);
	max_v = max_value(symbols, n);
	height = bits(max_v);

	OCC = new uint[max_v + 2];
	for (uint i = 0; i <= max_v + 1; i++)
		OCC[i] = 0;
	for (size_t i = 0; i < n; i++)
		OCC[symbols[i] + 1]++;

	size_t to_add = 0;
	for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0)
			to_add++;

	uint * new_symb = new uint[n + to_add];
	for (size_t i = 0; i < n; i++)
		new_symb[i] = symbols[i];

	if (deleteSymbols) {
		delete[] symbols;
		symbols = 0;
	}

	to_add = 0;
	for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0) {
			OCC[i]++;
			new_symb[n + to_add] = i - 1;
			to_add++;
		}

	size_t new_n = n + to_add;
	for (uint i = 1; i <= max_v + 1; i++)
		OCC[i] += OCC[i - 1];
	this->n = new_n;

	uint ** _bm = new uint*[height];
	for (uint i = 0; i < height; i++) {
		_bm[i] = new uint[new_n / W + 1];
		for (size_t j = 0; j < new_n / W + 1; j++)
			_bm[i][j] = 0;
	}

	build_level(_bm, new_symb, new_n, NULL);
	bitstring = new BitSequence*[height];
	C = new uint[height];
	for (uint i = 0; i < height; i++) {
		bitstring[i] = bmb->build(_bm[i], new_n);
		C[i] = bitstring[i]->rank0(new_n - 1);
		// cout << "C=" << C[i] << endl;
		delete[] _bm[i];
	}
	delete[] _bm;

	if (!deleteSymbols)
		for (size_t i = 0; i < n; i++)
			symbols[i] = am->unmap(symbols[i]);

	// delete [] new_symb; // already deleted in build_level()!
	// delete [] oc;
	bmb->unuse();
	// for(uint i=0;i<height;i++)
	//     cout << "C=" << C[i] << endl;
	this->length = n;
}

WaveletMatrix::WaveletMatrix() :
		Sequence(0) {
	bitstring = NULL;
	// occ = NULL;
	am = NULL;
}

WaveletMatrix::~WaveletMatrix() {
	if (bitstring) {
		for (uint i = 0; i < height; i++)
			if (bitstring[i])
				delete bitstring[i];
		delete[] bitstring;
	}
	// if(occ)
	// 	delete occ;
	if (am)
		am->unuse();
	delete[] C;
}

void WaveletMatrix::save(ofstream & fp) const {
	uint wr = WVMATRIX_HDR;
	saveValue(fp, wr);
	saveValue<size_t>(fp, n);
	saveValue(fp, max_v);
	saveValue(fp, height);
	saveValue(fp, C, height);
	am->save(fp);
	for (uint i = 0; i < height; i++)
		bitstring[i]->save(fp);	// XXX 没有实现
	// occ->save(fp);
	saveValue<uint>(fp, OCC, max_v + 2);
}

WaveletMatrix * WaveletMatrix::load(ifstream & fp) {
	uint rd = loadValue<uint>(fp);
	if (rd != WVMATRIX_HDR)
		return NULL;
	WaveletMatrix * ret = new WaveletMatrix();
	ret->n = loadValue<size_t>(fp);
	ret->length = ret->n;
	ret->max_v = loadValue<uint>(fp);
	ret->height = loadValue<uint>(fp);
	ret->C = loadValue<uint>(fp, ret->height);
	ret->am = Mapper::load(fp);
	if (ret->am == NULL) {
		delete ret;
		return NULL;
	}
	ret->am->use();
	ret->bitstring = new BitSequence*[ret->height];
	for (uint i = 0; i < ret->height; i++)
		ret->bitstring[i] = NULL;
	for (uint i = 0; i < ret->height; i++) {
		ret->bitstring[i] = BitSequence::load(fp);
		if (ret->bitstring[i] == NULL) {
			cout << "damn" << i << " " << ret->height << endl;
			delete ret;
			return NULL;
		}
	}
	ret->OCC = loadValue<uint>(fp, ret->max_v + 2);
	// ret->occ = BitSequence::load(fp);
	// if(ret->occ==NULL) {
	// 	delete ret;
	// 	return NULL;
	// }
	return ret;
}

inline uint get_start(uint symbol, uint mask) {
	return symbol & mask;
}

inline uint get_end(uint symbol, uint mask) {
	return get_start(symbol, mask) + !mask + 1;
}

bool WaveletMatrix::is_set(uint val, uint ind) const {
	assert(ind < height);
	return (val & (1 << (height - ind - 1))) != 0;
}

uint WaveletMatrix::set(uint val, uint ind) const {
	assert(ind <= height);
	return val | (1 << (height - ind - 1));
}

static size_t *x_start_left;
static size_t *x_start_right;
static size_t *x_end_left;
static size_t *x_end_right;

vector<uint> WaveletMatrix::n_range_intersect_aux(size_t *x_start,
		size_t *x_end, size_t n_ranges) {
	vector<uint> result;
	x_start_left = new size_t[n_ranges * height];
	x_start_right = new size_t[n_ranges * height];
	x_end_left = new size_t[n_ranges * height];
	x_end_right = new size_t[n_ranges * height];

	n_range_intersect(0, x_start, x_end, 0, n_ranges, x_start_left,
			x_start_right, x_end_left, x_end_right, result);

	delete x_start_left; // OJO se escribe asi?
	delete x_start_right; // OJO se escribe asi?
	delete x_end_left; // OJO se escribe asi?
	delete x_end_right; // OJO se escribe asi?
	return result;
}

void WaveletMatrix::n_range_intersect(uint lev, size_t * x_start, size_t *x_end,
		uint sym, size_t n_ranges, size_t *x_start_left, size_t *x_start_right,
		size_t *x_end_left, size_t *x_end_right, vector<uint> result) {

	// creo que estos tests estan de mas
	//	// llegamos a una hoja vacia
	//	if (start > end)
	//		return;
	//	// nos fuimos a la cresta
	//	if (end > n-1)
	//		return;
	//
	//	// revisamos que ningun rango se haya hecho vacio
	//	for (size_t i = 0 ; i < n_ranges;i++)
	//	{
	//		if (x_start[i] >  x_end[i])
	//			return;
	//	}

	// si aun no hemos llegado a la hoja
	//删去多余注释，除去无效传入参数，加入左右子树数量标识
	if (lev < height) {

		// construimos los arreglos para mantener los rangos

		size_t start_new = 0;
		size_t end_new;

		int liveleft = 1, liveright = 1;

		BitSequence* bs = bitstring[lev];

		// calculamos para cada rango, su nuevo intervalo
		for (size_t i = 0; i < n_ranges; i++) {
			start_new = x_start[i];
			if (x_start[i] > 0) {
				start_new = bs->rank1(start_new - 1);
				x_start_left[i] = x_start[i] - start_new /*+ 1*/;
				//XXX 有必要+1吗？
			} else
				x_start_left[i] = 0;

			end_new = bs->rank1(x_end[i]);
			if (liveleft) {
				x_end_left[i] = x_end[i] - end_new /*+ 1*/;
//				cout << "start_left:" << i << " " << x_start_left[i] << endl;
//				cout << "end_left:" << i << " " << x_end_left[i] << endl;
				//XXX 有必要+1吗？
				if (x_end_left[i] < x_start_left[i] || x_end[i] < end_new) {
					if (!liveright)
						//现在已经没有左子树了，右子树如果也没有
						return;
					liveleft = 0;
				}
			}
			if (liveright) {
				x_start_right[i] = start_new + C[lev];
				x_end_right[i] = end_new + C[lev] - 1;
//				cout << "start_right:" << i << " " << x_start_right[i] << endl;
//				cout << "end_right:" << i << " " << x_end_right[i] << endl;
				//end_new和C都是个数，表示序列得-1

				if (x_end_right[i] < x_start_right[i]) {
					if (!liveleft)
						return;
					liveright = 0;
				}
			}

		}
		// estas 2 no las optimizo porque en el wmatrix no estan

		// si no nos pasamos, ejecutamos a la izquierda
		if (liveleft) {
			n_range_intersect(lev + 1, x_start_left, x_end_left, sym, n_ranges,
					x_start_left + n_ranges, x_start_right + n_ranges,
					x_end_left + n_ranges, x_end_right + n_ranges, result);
//				start_left = start;
//				end_left = C[lev] - 1;
		}
		// si no nos pasamos, ejecutamos a la derecha

		if (liveright) {
			sym = sym | (1 << lev);
//				end_right = n - 1;
//				start_right = C[lev];
			n_range_intersect(lev + 1, x_start_right, x_end_right, sym,
					n_ranges, x_start_left + n_ranges, x_start_right + n_ranges,
					x_end_left + n_ranges, x_end_right + n_ranges, result);
		}
	}
	// llegamos a una hoja
	else {
		cout << "Adding sym = " << sym << endl;
		result.push_back(sym);
	}
}

vector<uint> WaveletMatrix::range_report_aux(size_t x_start, size_t x_end) {
	vector<uint> result;
//		result.reserve(1500);
	size_t num = x_end - x_start + 1;
	range_report(0, x_start, x_end, 0, num, result);
	return result;
}
void WaveletMatrix::range_report(uint lev, size_t x_start, size_t x_end,
		uint sym, size_t num, vector<uint> &result) {
	if (num <= 0)
		return;
	size_t x_start_left = 0, x_start_right = 0, x_end_left = 0, x_end_right = 0;

	size_t start_new = 0;
	size_t end_new = 0;

	if (lev < height) {
		BitSequence* bs = bitstring[lev];
		//x_start_left和x_end_left表示该层左子树的区域
		start_new = x_start;
		if (x_start > 0) {
			start_new = bs->rank1(start_new - 1);
			x_start_left = x_start - start_new;
		}

		end_new = bs->rank1(x_end);
		size_t num_right = (end_new - start_new);
		size_t num_left = num - num_right;

		x_end_left = x_end - end_new;
		// x_end == end_new ? x_end : x_end - end_new + 1;

		if (x_start_left <= x_end_left)
			range_report(lev + 1, x_start_left, x_end_left, sym, num_left,
					result);
//		sym = sym | (1 << (height - lev - 1));
		//进入右子树，开始赋值
		sym = sym | (1 << lev);

		x_start_right = start_new + C[lev];
		x_end_right = end_new + C[lev] - 1;

		if (x_start_right <= x_end_right)
			range_report(lev + 1, x_start_right, x_end_right, sym, num_right,
					result);
	} else {
		if (x_start > x_end)
			return;
		result.push_back(sym);
		cout << "Adding -> " << sym << " , " << endl;
//		result.push_back(am->unmap(sym));
//		cout << "Adding -> " << am->unmap(sym) << " , " << endl;

	}
}

uint WaveletMatrix::access(size_t pos) const {
	uint ret = 0;
	for (uint level = 0; level < height; level++) {
		size_t optR = 0;
		if (bitstring[level]->access(pos, optR)) {
			//进入右子树，可以赋值
			pos = C[level] + optR - 1;
			ret = ret | (1 << level);
		} else {
			//左子树不赋值
			pos = optR - 1;
		}
	}
	return ret;
//	return am->unmap(ret);
}

//表示符号symbol到pos位置（起点是0）所发生的次数
size_t WaveletMatrix::rank(uint symbol, size_t pos) const {
//	symbol = am->map(symbol);
	size_t start = 0;
	for (uint level = 0; level < height; level++) {
		if (is_set(symbol, height - level - 1)) {
			//进入右子树
			if (start > 0)
				start = bitstring[level]->rank1(start - 1);
			start += C[level];
			pos = bitstring[level]->rank1(pos) + C[level] - 1;
		} else {
			//进入左子树
			if (start > 0)		//这一步不同于waveletTree，因为1全部右移了
				start = start - bitstring[level]->rank1(start - 1);
			pos = pos - bitstring[level]->rank1(pos);
		}
		if (pos + 1 - start == 0)
			return 0;
	}
	return pos + 1 - start;
}

size_t WaveletMatrix::select(uint symbol, size_t j) const {
//	symbol = am->map(symbol);
	size_t pos = OCC[symbol] + j - 1; //(symbol == 0? -1 : occ->select1(symbol)) + j;
	for (int level = height - 1; level >= 0; level--) {
		// left
		if (pos < C[level]) {
			pos = bitstring[level]->select0(pos + 1);
		}					 // right
		else {
			pos = bitstring[level]->select1(pos - C[level] + 1);
		}
	}
	return pos;
}

size_t WaveletMatrix::getSize() const {
	size_t ptrs = sizeof(WaveletMatrix) + height * sizeof(Sequence*);
	size_t bytesBitstrings = 0;
	for (uint i = 0; i < height; i++)
		bytesBitstrings += bitstring[i]->getSize();
	return bytesBitstrings /* + occ->getSize() */+ ptrs + height * sizeof(uint)
			+ sizeof(uint) * (max_v + 2);
}

void WaveletMatrix::build_level(uint **bm, uint *symbols, size_t length,
		uint *occs) {
	uint sigma = max_value(symbols, length);					//记录length长度内最大值
	uint *new_order = new uint[sigma + 1];
	for (uint level = 0; level < height; level++) {
		/*for (uint i = 0 ; i < n;i++)
		 {
		 cout << " " << symbols[i];
		 }*/
		cout << endl;
		size_t zeroes = 0;
		for (uint i = 0; i < sigma + 1; i++)
			if (!is_set(i, height - level - 1)) {
				new_order[i] = 0;
			} else {
				new_order[i] = 1;
			}
		for (size_t i = 0; i < length; i++)
			if (!new_order[symbols[i]])
				zeroes++;
		uint *new_symbols = new uint[length];
		size_t new_pos0 = 0, new_pos1 = zeroes;
		for (size_t i = 0; i < length; i++) {
			if (!new_order[symbols[i]]) {
				new_symbols[new_pos0++] = symbols[i];
				bitclean(bm[level], i);
			} else {
				new_symbols[new_pos1++] = symbols[i];
				bitset(bm[level], i);
			}
		}
		delete[] symbols;
		symbols = new_symbols;
	}
	/*for (uint i = 0 ; i < n;i++)
	 {
	 cout << " " << symbols[i];
	 }
	 cout << endl; */

	delete[] symbols;
	delete[] new_order;
}

uint WaveletMatrix::max_value(uint *symbols, size_t n) {
	uint max_v = 0;
	for (size_t i = 0; i < n; i++)
		max_v = max(symbols[i], max_v);
	return max_v;
}

uint WaveletMatrix::bits(uint val) {
	uint ret = 0;
	while (val != 0) {
		ret++;
		val >>= 1;
	}
	return ret;
}

}
;
