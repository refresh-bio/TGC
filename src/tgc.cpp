#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <thread>
#include <list>
#include <utility>
#include <thread>

#include "port.h"
#include "qsmodel.h"
#include "rangecod.h"

#define MIN_MATCH_LEN		5
#define MAX_MATCH_LEN_EXP	11
#define MAX_MATCH_LEN		((1 << MAX_MATCH_LEN_EXP) + MIN_MATCH_LEN - 1)

#define HASH_LEN1			11
#define HASH_LEN2			2
#define HASH_STEP1			4
#define HASH_STEP2			(MIN_MATCH_LEN - HASH_LEN2 + 1)
#define HT_FACTOR			1.5

#define RC_SHIFT_REF			15
#define RC_SHIFT_FLAGS			15
#define RC_SHIFT_LITERALS		15
#define RC_SHIFT_LOG_LENS		15
#define RC_SHIFT_SUF_LENS		15
#define RC_SHIFT_PRE_OFFSETS	15
#define RC_SHIFT_SUF_OFFSETS	15

#define RC_SCALE_REF			(1u << 10)
#define RC_SCALE_FLAGS			(1u << 10)
#define RC_SCALE_LITERALS		(1u << 10)
#define RC_SCALE_LOG_LENS		(1u << 10)
#define RC_SCALE_SUF_LENS		(1u << 10)
#define RC_SCALE_PRE_OFFSETS	(1u << 10)
#define RC_SCALE_SUF_OFFSETS	(1u << 10)

/*#define RC_SCALE_REF			(1u << 13)
#define RC_SCALE_FLAGS			(1u << 7)
#define RC_SCALE_LITERALS		(1u << 13)
#define RC_SCALE_LOG_LENS		(1u << 9)
#define RC_SCALE_SUF_LENS		(1u << 11)
#define RC_SCALE_PRE_OFFSETS	(1u << 9)
#define RC_SCALE_SUF_OFFSETS	(1u << 11)
*/

#define RC_CTX_REF				4						// length of context for reference seq. compression

#define RC_CTX_LEN_FLAGS		9						// length of context for flags
#define RC_CTX_LEN_LITERALS		8						// length of context for literals
#define RC_MASK_LEN_LITERALS	0xFFFF					// mask for literals context

#define RC_CTX_BUCKET_LENS		(MAX_MATCH_LEN_EXP+1)	// no. of. contexts for prefixes of lengths

#define RC_CACHE_OFFSETS		4						// no. of cached offsets

//#define TURN_OFF_REF
//#define TURN_OFF_FLAGS
//#define TURN_OFF_LITERALS
//#define TURN_OFF_LENS
//#define TURN_OFF_OFFSETS

using namespace std;

//typedef int file_id_t;  //up to INT_MAX input BV files
typedef short file_id_t; //up to SHRT_MAX input BV files

bool compress_mode;

FILE *rc_file;			// output file (range coder compressed)
FILE *desc_file;		// output file (description of compressed data)
FILE *out_file;			// decompressed file

file_id_t *ht1, *ht2;
file_id_t *ht_zeros1, *ht_zeros2;

unsigned int no_files;							// no. of files to compress
unsigned long long ht_size1, ht_size2;			// size of hash table
unsigned long ht_slots1, ht_slots2;				// no. of slots in HT
unsigned long ht_slot_size_exp;					// size of hast table for each location
unsigned long ht_slot_size_mask;				// size of hast table for each location
unsigned long ht_slot_size;						// size of slot of hash table
int file_size;							// size of file

vector<unsigned char*> data;
vector<string> file_names;
set<string> file_names_decompress;

// Range coder and models
rangecoder rc;
qsmodel 
	qsm_ref_seq[RC_CTX_REF*8+1],
	qsm_flags[1u << RC_CTX_LEN_FLAGS], 
	qsm_literals[RC_CTX_LEN_LITERALS*8+1],
	qsm_offsets_pref, 
	*qsm_offsets_suf,
	qsm_lens_log[RC_CTX_BUCKET_LENS],
	qsm_lens_suf[RC_CTX_BUCKET_LENS];
	
// Contexts
unsigned int ctx_ref, ctx_flags, ctx_literals, ctx_lens, ctx_offsets;

int pos;
unsigned char ref_literal;
int cur_id_file;
vector<int> last_offsets(RC_CACHE_OFFSETS);


// ***************************************************************
bool check_data(int argc, char *argv[]);
bool prepare_files(string output_name);
void close_files(void);

void prepare_ht(void);
unsigned char *read_file(string &name);
void compress(void);
void parse_file(unsigned char *d);
void store_ref_file(unsigned char *d);

void uncompress_file(unsigned char *d);
void unstore_ref_file(unsigned char *d);
void decompress(void);

void init_rc_models(int mode);
void delete_rc_models(void);

void store_literal(unsigned char c);
void store_match(int id_file, int len);
void store_start(string &name);
void decode_literal(unsigned char *d, int &i);
void decode_match(unsigned char *d, int &i);

inline unsigned long hash_fun(unsigned char *p, int pos, int ver);
void insert_into_ht(file_id_t file_id, unsigned char *p, int ver);
inline pair<int, int> find_match(unsigned char *p, int pos, int ver);
inline int pop_count(unsigned int x);

// ***************************************************************
// Range coder-relared functions
// ***************************************************************

// ***************************************************************
void save_byte(int x)
{
	putc(x, rc_file);
}

// ***************************************************************
int read_byte()
{
	int a = getc(rc_file);
	return a;
}


// ***************************************************************
// Initialize range coder models
void init_rc_models(int mode)
{
	// Reference sequence
	for(int i = 0; i < RC_CTX_REF*8+1; ++i)
		initqsmodel(&qsm_ref_seq[i], 256, RC_SHIFT_REF, RC_SCALE_REF, NULL, mode);

	// Flags
	for(int i = 0; i < (1u << RC_CTX_LEN_FLAGS); ++i)
		initqsmodel(&qsm_flags[i], 2, RC_SHIFT_FLAGS, RC_SCALE_FLAGS, NULL, mode);

	// Literals
	for(int i = 0; i < RC_CTX_LEN_LITERALS*8+1; ++i)
		initqsmodel(&qsm_literals[i], 257, RC_SHIFT_LITERALS, RC_SCALE_LITERALS, NULL, mode);
	
	// Lenghts
	for(int i = 0; i < RC_CTX_BUCKET_LENS; ++i)
		initqsmodel(&qsm_lens_log[i], RC_CTX_BUCKET_LENS, RC_SHIFT_LOG_LENS, RC_SCALE_LOG_LENS, NULL, mode);

	for(int i = 2; i < RC_CTX_BUCKET_LENS; ++i)
	{
		int size = 1u << (i-1);
		initqsmodel(&qsm_lens_suf[i], size, RC_SHIFT_SUF_LENS, RC_SCALE_SUF_LENS, NULL, mode);
	}
		
	// Offsets
	int rc_ctx_bucket_offsets = (no_files + 255) / 256 + 1;
	qsm_offsets_suf = new qsmodel[rc_ctx_bucket_offsets];
	
	initqsmodel(&qsm_offsets_suf[0], RC_CACHE_OFFSETS, RC_SHIFT_SUF_OFFSETS, RC_SCALE_SUF_OFFSETS, NULL, mode);

	for(int i = 0; i < rc_ctx_bucket_offsets-1; ++i)
	{
		int size = 256;
		if(i * 256 + 256 > no_files)
			size = (no_files) - i * 256;
		initqsmodel(&qsm_offsets_suf[i+1], size, RC_SHIFT_SUF_OFFSETS, RC_SCALE_SUF_OFFSETS, NULL, mode);
	}
	initqsmodel(&qsm_offsets_pref, rc_ctx_bucket_offsets, RC_SHIFT_PRE_OFFSETS,	RC_SCALE_PRE_OFFSETS, NULL, mode);
}

// ***************************************************************
// Delete range coder models
void delete_rc_models(void)
{
	// Reference sequence
	for(int i = 0; i < RC_CTX_REF*8+1; ++i)
		deleteqsmodel(&qsm_ref_seq[i]);

	// Flags
	for(int i = 0; i < (1u << RC_CTX_LEN_FLAGS); ++i)
		deleteqsmodel(&qsm_flags[i]);

	// Literals
	for(int i = 0; i < RC_CTX_LEN_LITERALS*8+1; ++i)
		deleteqsmodel(&qsm_literals[i]);
	
	// Lenghts
	for(int i = 0; i < RC_CTX_BUCKET_LENS; ++i)
		deleteqsmodel(&qsm_lens_log[i]);

	for(int i = 2; i < RC_CTX_BUCKET_LENS; ++i)
		deleteqsmodel(&qsm_lens_suf[i]);
		
	// Offsets
	int rc_ctx_bucket_offsets = (no_files + 255) / 256 + 1;
	for(int i = 0; i < rc_ctx_bucket_offsets; ++i)
		deleteqsmodel(&qsm_offsets_suf[i]);

	delete[] qsm_offsets_suf;
	deleteqsmodel(&qsm_offsets_pref);
}


// ***************************************************************
// Utility functions
// ***************************************************************

// **************************************************************
// Calculate number of set bits
inline int pop_count(unsigned int x)
{
	int cnt = 0;
	
	for(; x; ++cnt)
		x &= x-1;
		
	return cnt;
}


// ***************************************************************
// Hash table functions
// ***************************************************************

// **************************************************************
inline unsigned long hash_fun(unsigned char *p, int pos, int ver)
{
	unsigned long res = 0;
	
	if(ver == 1)
		for(int i = 0; i < HASH_LEN1; ++i)
			res = res * 65531u + p[pos+i] * 29u;
	else
		for(int i = 0; i < HASH_LEN2; ++i)
			res = res * 65531u + p[pos+i] * 29u;
		
	return res & ht_slot_size_mask;
}

// ***************************************************************
inline pair<int, int> find_match(unsigned char *p, int pos, int ver)
{
	int best_id  = -1;
	int best_len = 0;
	int hash_step, hash_len;
	
	if(ver == 1)
	{
		hash_step = HASH_STEP1;
		hash_len  = HASH_LEN1;
	}
	else
	{
		hash_step = HASH_STEP2;
		hash_len  = HASH_LEN2;
	}
			
	int pos_norm = (pos + hash_step - 1) / hash_step * hash_step;
	
	if(pos_norm + hash_len > file_size)
		return make_pair(-1, -1);
		
	unsigned long h = hash_fun(p, pos_norm, ver);
	unsigned long long off = ((unsigned long long)pos_norm / hash_step) << (ht_slot_size_exp);
	file_id_t *ht;
	
	if(ver == 1)
		ht = ht1;
	else
		ht = ht2;
	
	while(ht[off+h] != -1)
	{
		int file_id = ht[off+h];
		int i;
		
		if(file_id < cur_id_file)
		{
			for(i = pos; i < file_size; ++i)
			{
				if(p[i] != data[file_id][i])
					break;
				if(++i >= file_size)
					break;
				if(p[i] != data[file_id][i])
					break;
			}

			if(i - pos > best_len)
			{
				best_len = i - pos;
				best_id  = file_id;
			}	
		}
		
		h = (h+1) & ht_slot_size_mask;
	}
	
	if(best_len < MIN_MATCH_LEN)
		return make_pair(-1, -1);
		
	return make_pair(best_id, best_len);
}

// ***************************************************************
void insert_into_ht(file_id_t file_id, unsigned char *p, int ver)
{
	unsigned int n_slot = 0;
	unsigned long long off = 0;
	
	int hash_len, hash_step;
	file_id_t *ht, *ht_zeros;
	
	if(ver == 1)
	{
		hash_len  = HASH_LEN1;
		hash_step = HASH_STEP1;
		ht        = ht1;
		ht_zeros  = ht_zeros1;
	}
	else
	{
		hash_len  = HASH_LEN2;
		hash_step = HASH_STEP2;
		ht        = ht2;
		ht_zeros  = ht_zeros2;
	}
	
	for(int i = 0; i < file_size - hash_len + 1; i += hash_step, ++n_slot, off += ht_slot_size)
	{
		unsigned long h = hash_fun(p, i, ver);

		if(h <= ht_zeros[n_slot])		// insert counter of entries with hash_value = 0
			h = ht_zeros[n_slot]++;
	
		while(ht[off+h] != -1)
			h = (h+1) & ht_slot_size_mask;

		ht[off+h] = file_id;
	}
}

// ***************************************************************
void prepare_ht(void)
{
	for(ht_slot_size_exp = 5; no_files * HT_FACTOR > (1u << ht_slot_size_exp); ht_slot_size_exp++)
		;
	ht_slot_size      = 1u << ht_slot_size_exp;
	ht_slot_size_mask = ht_slot_size-1;
	ht_slots1	      = (unsigned long) (file_size / HASH_STEP1 + 1);
	ht_slots2	      = (unsigned long) (file_size / HASH_STEP2 + 1);
	ht_size1          = ht_slot_size * ht_slots1;
	ht_size2          = ht_slot_size * ht_slots2;
	
	cout << "HT sizes: " << ht_size1 / (1<<20)*sizeof(file_id_t) + ht_size2 / (1<<20)*sizeof(file_id_t) << "MB\n";

	ht1 = new file_id_t[ht_size1];
	fill(ht1, ht1+ht_size1, -1);
	ht2 = new file_id_t[ht_size2];
	fill(ht2, ht2+ht_size2, -1);

	// set counter of entries with hash value = 0
	ht_zeros1 = new file_id_t[ht_slots1];
	fill(ht_zeros1, ht_zeros1+ht_slots1, 0);
	ht_zeros2 = new file_id_t[ht_slots2];
	fill(ht_zeros2, ht_zeros2+ht_slots2, 0);
}


// ***************************************************************
// I/Ofunctions
// ***************************************************************

// ***************************************************************
unsigned char* read_file(string &name)
{
	FILE *in = fopen(name.c_str(), "rb");
	if(!in)
	{
		cout << "No file: " << name << "\n";
		return NULL;
	}

	// Check size
	fseek(in, 0, SEEK_END);
	if(ftell(in) != (unsigned int) file_size)
	{
		cout << "File " << name << " is of incompatibile size\n";
		fclose(in);
		return NULL;
	}
	fseek(in, 0, SEEK_SET);

	unsigned char *d = new unsigned char[file_size];
	fread(d, 1, file_size, in);
	fclose(in);

	return d;
}

// ***************************************************************
// Prepare output files
bool prepare_files(string output_name)
{
	if(compress_mode)
	{
		rc_file   = fopen((output_name + ".tgc_data").c_str(), "wb");
		desc_file = fopen((output_name + ".tgc_desc").c_str(), "wt");
	
		fprintf(desc_file, "%d\n", file_size);
	
		return rc_file && desc_file;
	}
	else
	{
		rc_file   = fopen((output_name + ".tgc_data").c_str(), "rb");
		
		return rc_file;
	}
}

// ***************************************************************
// Close output files
void close_files(void)
{
	if(rc_file)
		fclose(rc_file);
	if(desc_file)
		fclose(desc_file);
	if(out_file)
		fclose(out_file);
}


// ***************************************************************
// Compression functions
// ***************************************************************

// ***************************************************************
// Store description of the processed file
void store_start(string &name)
{
	fprintf(desc_file, "%s\n", name.c_str());
	cout << "\rProcessing " << name << " (" << cur_id_file+1 << " of " << file_names.size() << ")    ";
}

// ***************************************************************
void store_literal(unsigned char c)
{
	int syfreq, ltfreq;
	
	// Flag
	qsgetfreq(&qsm_flags[ctx_flags], 0, &syfreq, &ltfreq);
#ifndef TURN_OFF_FLAGS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_FLAGS);
#endif
	qsupdate(&qsm_flags[ctx_flags], 0);
	ctx_flags = ((ctx_flags << 1) + 0) & ((1u << RC_CTX_LEN_FLAGS) - 1);

	int ci = c;
	if(ci == 0)
		ci = 256;
	else
		ci ^= ref_literal;
	
	int cnt = pop_count(ctx_literals) + 16 * (ctx_flags & 0x3);

	qsgetfreq(&qsm_literals[cnt], ci, &syfreq, &ltfreq);
#ifndef TURN_OFF_LITERALS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LITERALS);
#endif
	qsupdate(&qsm_literals[cnt], ci);
	
	ctx_literals = ((ctx_literals << 8) + c) & RC_MASK_LEN_LITERALS;	
	ref_literal  = 0;
}

// ***************************************************************
void store_match(int id_file, int len)
{
	int syfreq, ltfreq;

	// Flags
	qsgetfreq(&qsm_flags[ctx_flags], 1, &syfreq, &ltfreq);
#ifndef TURN_OFF_FLAGS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_FLAGS);
#endif
	qsupdate(&qsm_flags[ctx_flags], 1);
	ctx_flags = ((ctx_flags << 1) + 1) & ((1u << RC_CTX_LEN_FLAGS) - 1);
	
	if(pos + len < file_size)
		ref_literal = data[id_file][pos+len];
	else
		ref_literal = 0;

	// Lens
	len -= MIN_MATCH_LEN;
	int p1, p2;
	
	if(len < 2)
		p1 = len;
	else 
	{
		for(p1 = 2; len >= (1u << p1); ++p1)
			;
		p2 = len - (1u << (p1-1));
	}
	
	qsgetfreq(&qsm_lens_log[ctx_lens], p1, &syfreq, &ltfreq);
#ifndef TURN_OFF_LENS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_LOG_LENS);
#endif
	qsupdate(&qsm_lens_log[ctx_lens], p1);
	
	if(p1 > 1)
	{
		qsgetfreq(&qsm_lens_suf[p1], p2, &syfreq, &ltfreq);
#ifndef TURN_OFF_LENS
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SUF_LENS);
#endif
		qsupdate(&qsm_lens_suf[p1], p2);
	}
	ctx_lens = ctx_lens * 0.75 + p1 * 0.25;
	if(ctx_lens >= RC_CTX_BUCKET_LENS)
		ctx_lens = RC_CTX_BUCKET_LENS-1;

	// ***** Offsets
	int p;
	for(p = 0; p < RC_CACHE_OFFSETS; ++p)
		if(last_offsets[p] == id_file)
			break;
			
	if(p < RC_CACHE_OFFSETS)
	{
		if(p)
		{
			last_offsets.erase(last_offsets.begin()+p);
			last_offsets.insert(last_offsets.begin(), id_file);
		}
	}
	else
	{
		last_offsets.insert(last_offsets.begin(), id_file);
		last_offsets.resize(RC_CACHE_OFFSETS);
	}
	
	int of1, of2;
	if(p < RC_CACHE_OFFSETS)
	{
		of1 = 0;
		of2 = p;
	}
	else
	{
		of1 = id_file / 256 + 1;
		of2 = id_file % 256;
	}
	
	qsgetfreq(&qsm_offsets_pref, of1, &syfreq, &ltfreq);
#ifndef TURN_OFF_OFFSETS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_PRE_OFFSETS);
#endif
	qsupdate(&qsm_offsets_pref, of1);
	
	qsgetfreq(&qsm_offsets_suf[of1], of2, &syfreq, &ltfreq);
#ifndef TURN_OFF_OFFSETS
	encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_SUF_OFFSETS);
#endif
	qsupdate(&qsm_offsets_suf[of1], of2);	
}

// ***************************************************************
void parse_file(unsigned char * d)
{
	int n_files = data.size();
	int best_len;
	int best_id;
	
	pos = 0;
	while(pos < file_size)
	{
		pair<int, int> match = find_match(d, pos, 1);
		if(match.second < HASH_LEN1 + HASH_STEP1 - 1)
			match = find_match(d, pos, 2);
		
		if(match.first < 0)
		{
			store_literal(d[pos]);
			pos++;
		}
		else
		{
			if(match.second > MAX_MATCH_LEN)
				match.second = MAX_MATCH_LEN;
			store_match(match.first, match.second);
			pos += match.second;
		}
	}
}

// ***************************************************************
// Compress reference file
void store_ref_file(unsigned char *d)
{
#ifdef TURN_OFF_REF
	return;
#endif

	int syfreq, ltfreq;

	for(int i = 0; i < file_size; ++i)
	{
		// Calculate no. of bits in ctx
		int c = pop_count(ctx_ref);

		qsgetfreq(&qsm_ref_seq[c], d[i], &syfreq, &ltfreq);
		encode_shift(&rc, syfreq, ltfreq, RC_SHIFT_REF);
		qsupdate(&qsm_ref_seq[c], d[i]);
		
		ctx_ref = (ctx_ref << 8) + d[i];
	}
}

// ***************************************************************
void compress(void)
{
	unsigned char *d;
	
	init_rc_models(1);					// compression mode (1)
	start_encoding(&rc, 0, 0);
	last_offsets.resize(RC_CACHE_OFFSETS);

	for(int i = 0; i < file_names.size(); ++i)
	{
		if((d = read_file(file_names[i])) == NULL)
			continue;
		
		store_start(file_names[i]);

		data.push_back(d);
		thread t1([&]{insert_into_ht(cur_id_file, d, 1);});
		thread t2([&]{insert_into_ht(cur_id_file, d, 2);});

		if(i)
			parse_file(d);
		else
			store_ref_file(d);
			
		t1.join();
		t2.join();
		++cur_id_file;
	}
	
	done_encoding(&rc);
	delete_rc_models();
}

// ***************************************************************
// Decompression functions
// ***************************************************************

// ***************************************************************
void decode_literal(unsigned char *d, int &i)
{
	int syfreq, ltfreq;
	int cnt = pop_count(ctx_literals) + 16 * (ctx_flags & 0x3);

	ltfreq = decode_culshift(&rc, RC_SHIFT_LITERALS);
	int ch = qsgetsym(&qsm_literals[cnt], ltfreq);

	unsigned char c = ch;
	
	if(ch == 256)
		c = 0;
	else
		c ^= ref_literal;
		
	qsgetfreq(&qsm_literals[cnt], ch, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_LITERALS);
	qsupdate(&qsm_literals[cnt], ch);

	d[i++] = c;
	
	ctx_literals = ((ctx_literals << 8) + c) & RC_MASK_LEN_LITERALS;
	ref_literal  = 0;
}

// ***************************************************************
void decode_match(unsigned char *d, int &i)
{
	int syfreq, ltfreq;
	int p1, p2;
	int ch;
	
	int len, id_file;
	int of1, of2, p;
	
	// ***** Lens
	ltfreq = decode_culshift(&rc, RC_SHIFT_LOG_LENS);
	p1 = qsgetsym(&qsm_lens_log[ctx_lens], ltfreq);
	qsgetfreq(&qsm_lens_log[ctx_lens], p1, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_LOG_LENS);
	qsupdate(&qsm_lens_log[ctx_lens], p1);
	
	if(p1 > 1)
	{
		ltfreq = decode_culshift(&rc, RC_SHIFT_SUF_LENS);
		p2 = qsgetsym(&qsm_lens_suf[p1], ltfreq);
		qsgetfreq(&qsm_lens_suf[p1], p2, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SUF_LENS);
		qsupdate(&qsm_lens_suf[p1], p2);
		
		len = p2 + (1u << (p1-1));
	}
	else
		len = p1;
		
	len += MIN_MATCH_LEN;

	ctx_lens = ctx_lens * 0.75 + p1 * 0.25;
	if(ctx_lens >= RC_CTX_BUCKET_LENS)
		ctx_lens = RC_CTX_BUCKET_LENS-1;

	// ***** Offsets
	ltfreq = decode_culshift(&rc, RC_SHIFT_PRE_OFFSETS);
	of1 = qsgetsym(&qsm_offsets_pref, ltfreq);
	qsgetfreq(&qsm_offsets_pref, of1, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_PRE_OFFSETS);
	qsupdate(&qsm_offsets_pref, of1);
	
	ltfreq = decode_culshift(&rc, RC_SHIFT_SUF_OFFSETS);
	of2 = qsgetsym(&qsm_offsets_suf[of1], ltfreq);
	qsgetfreq(&qsm_offsets_suf[of1], of2, &syfreq, &ltfreq);
	decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_SUF_OFFSETS);
	qsupdate(&qsm_offsets_suf[of1], of2);
	
	if(of1 == 0)
	{
		id_file = last_offsets[of2];
		p = of2;
	}
	else
	{
		id_file = (of1 - 1) * 256 + of2;
		p = RC_CACHE_OFFSETS;
	}
	
	if(p < RC_CACHE_OFFSETS)
	{
		if(p)
		{
			last_offsets.erase(last_offsets.begin()+p);
			last_offsets.insert(last_offsets.begin(), id_file);
		}
	}
	else
	{
		last_offsets.insert(last_offsets.begin(), id_file);
		last_offsets.resize(RC_CACHE_OFFSETS);
	}
	
	for(int j = 0; j < len; ++j, ++i)
		d[i] = data[id_file][i];

	if(i < file_size)
		ref_literal = data[id_file][i];
	else
		ref_literal = 0;
}
	
// ***************************************************************
void uncompress_file(unsigned char *d)
{
	int syfreq, ltfreq;
	int flag;
	
	for(int i = 0; i < file_size; )
	{
		// Flag
		ltfreq = decode_culshift(&rc, RC_SHIFT_FLAGS);
		flag = qsgetsym(&qsm_flags[ctx_flags], ltfreq);
		qsgetfreq(&qsm_flags[ctx_flags], flag, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_FLAGS);
		qsupdate(&qsm_flags[ctx_flags], flag);
		
		ctx_flags = ((ctx_flags << 1) + flag) & ((1u << RC_CTX_LEN_FLAGS) - 1);

		if(flag == 0)			// literal
			decode_literal(d, i);
		else if(flag == 1)		// match
			decode_match(d, i);			
		else
		{
			cout << "Corrupted file\n";
			exit(1);
		}
	}
}

// ***************************************************************
void unstore_ref_file(unsigned char *d)
{
	int ch;
	int syfreq, ltfreq;

	for(int i = 0; i < file_size; ++i)
	{
		// Calculate no. of bits in ctx
		int c = pop_count(ctx_ref);
			
		ltfreq = decode_culshift(&rc, RC_SHIFT_REF);
		ch = qsgetsym(&qsm_ref_seq[c], ltfreq);
		
		d[i] = ch;
		qsgetfreq(&qsm_ref_seq[c], ch, &syfreq, &ltfreq);
		decode_update(&rc, syfreq, ltfreq, 1u << RC_SHIFT_REF);
		qsupdate(&qsm_ref_seq[c], ch);
		
		ctx_ref = (ctx_ref << 8) + ch;
	}
}

// ***************************************************************
void decompress(void)
{
	init_rc_models(0);				// decompression mode (0)
	start_decoding(&rc);
	last_offsets.resize(RC_CACHE_OFFSETS);
	
	for(int i = 0; i < file_names.size(); ++i)
	{
		unsigned char *d = new unsigned char[file_size];
		
		if(i)
			uncompress_file(d);
		else
			unstore_ref_file(d);
		data.push_back(d);

		if(file_names_decompress.find(file_names[i]) != file_names_decompress.end())
		{
			string out_name = file_names[i] + ".ori";
			out_file = fopen(out_name.c_str(), "wb");
			if(!out_file)
			{
				cout << "Cannot create file " << out_name << "\n";
				exit(1);
			}
			
			fwrite(d, 1, file_size, out_file);
			fclose(out_file);
			
			file_names_decompress.erase(file_names[i]);
		}
		
		if(file_names_decompress.empty())
			break;
	}

	out_file = NULL;
	done_decoding(&rc);
}


// ***************************************************************
// Main program
// ***************************************************************

// ***************************************************************
// Check parameters, input files etc.
bool check_data(int argc, char *argv[])
{
	if((argc < 2) || (argc < 3 && strcmp(argv[1], "d") == 0) || (argc < 4 && strcmp(argv[1], "c") == 0))
	{
		cout << "Usage: tgc <mode> <output_name> [list_file_name]\n";
		cout << "  mode           - c (compression), d (decompression)\n";
		cout << "  output_name    - name of the output file\n";
		cout << "  list_file_name - name of the file with the list of files to compress or decompress\n";
		cout << "                   (in decompression if absent, all files will be decompressed)\n";
		cout << "Examples:\n";
		cout << "  tgc c chr1 files_chr1\n";
		cout << "  tgc d chr1\n";
		return false;
	}

	if(strcmp(argv[1], "c") == 0)		// compression
		compress_mode = true;
	else 
		compress_mode = false;

	if(compress_mode)
	{
		ifstream inf(argv[3]);
		istream_iterator<string> inf_iter(inf);
		file_names.assign(inf_iter, istream_iterator<string>());
	}
	else
	{
		ifstream inf(string(argv[2]) + ".tgc_desc");
		istream_iterator<string> inf_iter(inf);
		file_names.assign(inf_iter, istream_iterator<string>());
		
		file_size = atoi(file_names.front().c_str());
		file_names.erase(file_names.begin());

		if(argc > 3)
		{
			ifstream inf_dec(argv[3]);
			istream_iterator<string> inf_dec_iter(inf_dec);
			file_names_decompress.insert(inf_dec_iter, istream_iterator<string>());
		}
		else
			file_names_decompress.insert(file_names.begin(), file_names.end());
	}

	no_files = file_names.size();

	if(!no_files)
	{
		cout << "There are no files to process\n";
		return false;
	}

	if(compress_mode)		// compression
	{
		// Check size of the first file - it is assummed that all files have the same size!
		FILE *in = fopen(file_names[0].c_str(), "rb");
		if(!in)
		{
			cout << "The reference file " << file_names[0] << " does not exist\n";
			return false;
		}
	
		fseek(in, 0, SEEK_END);
		file_size = (int) ftell(in);
		fclose(in);
	}
	
	return true;
}

// ***************************************************************
// ***************************************************************
int main(int argc, char *argv[])
{
	clock_t t1 = clock();

	if(!check_data(argc, argv))
		return 0;
		
	if(!prepare_files(string(argv[2])))
		return 0;
	
	if(compress_mode)
	{
		cout << "Initializing\n";
		prepare_ht();
	
		cout << "Compressing...\n";
		compress();
	
		close_files();
	}
	else
	{
		cout << "Decompressing...\n";
		decompress();

		close_files();
	}
	
	cout << "\nTime: " << (double) (clock() - t1) / CLOCKS_PER_SEC << "\n";
	
	return 0;
}
