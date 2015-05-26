#include <cstdlib>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <math.h>

using namespace std;

const string MID_SEQ="GTGAGCGGATAACAAT";
const string END_SEQ="TGGAA";
const string QUAL_SCORE = " !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

bool tag_match(const string& dna_tag, const string& rna_tag) 
{
  if (rna_tag.size() > dna_tag.size())
  {
    if (rna_tag.substr(rna_tag.size() - dna_tag.size(), dna_tag.size()) == dna_tag)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  if (rna_tag == dna_tag.substr(dna_tag.size() - rna_tag.size(), rna_tag.size()))
  {
    return true;
  }
  else
  {
    return false;
  }
}


int main(int argc, char * argv[])
{
  if(argc < 10)
  {
    cerr << "Usage: " << argv[0] << " [quality score cutoff] [digital tag length] [tag length] [extra position] [key length] <output file for unformated result> <run stats file name> <dna parsed file> <rna FASTQ files>\n";
    return 0;
  }

  const int qual_min = atoi(argv[1]);
  const int digital_tag_len = atoi(argv[2]), tag_len = atoi(argv[3]), extra_pos = atoi(argv[4]), key_len = atoi(argv[5]);
  const int app_tag_len = tag_len + extra_pos;
  const int total_pos = tag_len + extra_pos + 1;
  char qual_min_char = QUAL_SCORE[qual_min];
  string rna_file_name = argv[6];
  
  string sample_id = rna_file_name.substr(0,5);
  ifstream fastq_file, dna_fq_file(argv[8]);
  ofstream unformated_result(argv[6]), run_log(argv[7]), tag_record;
  tag_record.open((rna_file_name.substr(0, 7) + "_tag_records.txt").c_str());

  unformated_result << "\t\t" << (sample_id + "_match") << '\t' << (sample_id + "_digital_match") << '\t' << (sample_id + "_all") << '\t' << (sample_id + "_digital_all") << '\t' << (sample_id + "_match_percent") << '\t' << (sample_id + "_digital_match_percent") << '\t' << (sample_id + "_all_percent") << '\t' <<(sample_id + "_all_percent") << '\t'<< (sample_id + "_mismatch_percent") << '\t' << (sample_id + "_digital_mismatch_percent") << '\n';
  
  string file_name;
  string label_hold,bases,pos,qual;
  int flag = 0;

  int mid_index = 0;
  int q_fail = 0, c_fail = 0;
  int total_reads = 0;
  int parsed_reads = 0;
  int discard_reads = 0;
  int recovered_reads = 0;

  string dna_key, dna_tag, rna_key, rna_tag, rna_dig_tag;
  int dna_key_count;

  set<string> raw_rna_digital_tag_set, raw_rna_key_set, raw_rna_tag_set;
  map<string, int> raw_rna_key_map, raw_rna_tag_map, raw_rna_digital_tag_map;
  map<string, pair<string, int> > dna_key_map;
  map<string, map<string, int> > dna_tag_map;
  map<string, map<string, pair<int, set<string> > > > rna_map;//{key:{tag:(count, dig_tag_set),...},...}
  map<string, set<string> > rna_tag_map;
  map<string, int> rna_key_count_map, rna_digital_tag_key_count_map;
  map<string, set<string> > rna_digital_tag_key_div_map;

  int no_dna_key_match = 0, good_reads_count = 0;
  set<string> no_dna_key_match_set, rna_good_key_set;


  while (dna_fq_file)
  {
    dna_fq_file >> dna_key >> dna_tag >> dna_key_count;
    if (dna_fq_file.eof())
    {
      break;
    }

    //cout << dna_key << ' ' << dna_tag << ' ' << dna_key_count << endl;
    dna_key_map[dna_key].first = dna_tag;
    dna_key_map[dna_key].second = dna_key_count;
    dna_tag_map[dna_tag][dna_key] = dna_key_count;
    if (dna_tag_map[dna_tag].find("SUM") == dna_tag_map[dna_tag].end()) 
    {
      dna_tag_map[dna_tag]["SUM"] = dna_key_count;
    }
  }
  dna_fq_file.close();
  
  for(int i = 9; i < argc; i++)
  {
    file_name = argv[i];
    fastq_file.open(file_name.c_str());
    if(!fastq_file)
    {
      cerr << "Error reading file\n";
      return 0;
    }
  
    while(fastq_file)
    {
      getline(fastq_file,label_hold);
      getline(fastq_file,bases);
      getline(fastq_file,pos);
      getline(fastq_file,qual);
      if(fastq_file.eof())
	break;
      total_reads++;
      flag = 0;
      mid_index = bases.find(MID_SEQ.substr(extra_pos, MID_SEQ.length() - extra_pos)) - extra_pos;

      if(mid_index < 0 || mid_index + 21 + key_len > bases.length() || mid_index - digital_tag_len + extra_pos < 1)
      {
	c_fail++; //Size fails.
	//q_fail++;
	flag = 1;
      }
      else
      {
	rna_key = bases.substr(mid_index + MID_SEQ.length(), key_len);
	//cout << rna_key << endl;
	rna_dig_tag = bases.substr(0, digital_tag_len);
	rna_tag = bases.substr(digital_tag_len, mid_index - digital_tag_len + extra_pos);
	raw_rna_tag_set.insert(rna_tag);
	raw_rna_key_set.insert(rna_key);
	raw_rna_digital_tag_set.insert(rna_dig_tag);
	recovered_reads++;
	for(int j = 0; j < (mid_index + key_len + 21); j++)
	{
	  if((qual.at(j) < qual_min_char) || bases.at(j) == 'N')
	  {
	    q_fail++; //Qual fails.
	    flag = 1;
	    break;
	  }
	}

	if (END_SEQ.compare(bases.substr(mid_index + key_len + 16, END_SEQ.length())) != 0)
	{
	  flag = 1;
	  //if (flag != 1)
	  //{
	  //  q_fail++;
	  //}
	}
      }

      if(flag != 1)
      {
	rna_map[rna_key][rna_tag].first++;
	rna_map[rna_key][rna_tag].second.insert(rna_dig_tag);
	if (dna_key_map.find(rna_key) != dna_key_map.end())
	{
	  //cout << rna_key << endl;
	  rna_tag_map[dna_key_map[rna_key].first].insert(rna_key);
	  rna_key_count_map[dna_key_map[rna_key].first]++;
	  rna_digital_tag_key_count_map[rna_dig_tag]++;
	  rna_digital_tag_key_div_map[rna_dig_tag].insert(rna_key);
	  rna_good_key_set.insert(rna_key);
	  good_reads_count++;
	}
	else
	{
	  no_dna_key_match++;
	  no_dna_key_match_set.insert(rna_key);
	}	  
	parsed_reads++;
      }
      else
      {
	discard_reads++;
      }
    }
    fastq_file.close();
  }

  int all_matched_reads = 0, all_matched_digital_reads = 0;
  
  for (map<string, set<string> >::iterator it = rna_tag_map.begin(); it != rna_tag_map.end(); it++)
  {
    int match[total_pos], digital_match[total_pos], all[total_pos], digital_all[total_pos], mis[total_pos], dig_mis[total_pos];
    //int match[], digital_match[], all[], digital_all[], mis[], dig_mis[];
    memset( match, 0, total_pos * sizeof(int) );
    memset( digital_match, 0, total_pos * sizeof(int) );
    memset( all, 0, total_pos * sizeof(int) );
    memset( digital_all, 0, total_pos * sizeof(int) );
    memset( mis, 0, total_pos * sizeof(int) );
    memset( dig_mis, 0, total_pos * sizeof(int) );

    double sum_match = 0, sum_dig_match = 0, sum_mis = 0, sum_dig_mis = 0, sum_all = 0, sum_dig_all = 0;
    for (set<string>::iterator itsec = (it->second).begin(); itsec != (it->second).end(); itsec++)
    {
      for (map<string, pair<int, set<string> > >::iterator it_third = (rna_map[*itsec]).begin(); it_third != (rna_map[*itsec]).end(); it_third++)
      {
        int ic_length = (it_third->first).length();
	tag_record << (it_third->first) << '\t' << (it->first) << '\t' << (12 - ic_length) << '\t' << (int)((it_third->second).first) << '\t' << (int)((it_third->second).second.size()) << '\t';
	if ((it_third->first).size() > app_tag_len)
	{
	  if (tag_match(it->first, it_third->first))
	  {
	    match[0] += (it_third->second).first;
	    digital_match[0] += (it_third->second).second.size();
	    all[0] += (it_third->second).first;
	    digital_all[0] += (it_third->second).second.size();
	    all_matched_reads += (it_third->second).first;
	    sum_match += (it_third->second).first;
	    sum_dig_match += (it_third->second).second.size();
	    sum_all += (it_third->second).first;
	    sum_dig_all += (it_third->second).second.size();
	    all_matched_digital_reads += (it_third->second).second.size();
	    tag_record << '1' << endl;
	  }
	  else
	  {
	    mis[0] += (it_third->second).first;
	    dig_mis[0] += (it_third->second).second.size();
	    all[0] += (it_third->second).first;
	    digital_all[0] += (it_third->second).second.size();
	    sum_mis += (it_third->second).first;
	    sum_dig_mis += (it_third->second).second.size();
	    sum_all += (it_third->second).first;
	    sum_dig_all += (it_third->second).second.size();
	    tag_record << '0' << endl;
	  }
	}
	else
	{
	  int pos = app_tag_len + 1 - (it_third->first).size();
	  if (tag_match(it->first, it_third->first))
	  {
	    match[pos] += (it_third->second).first;
	    digital_match[pos] += (it_third->second).second.size();
	    all[pos] += (it_third->second).first;
	    digital_all[pos] += (it_third->second).second.size();
	    all_matched_reads += (it_third->second).first;
	    sum_match += (it_third->second).first;
	    sum_dig_match += (it_third->second).second.size();
	    sum_all += (it_third->second).first;
	    sum_dig_all += (it_third->second).second.size();
	    all_matched_digital_reads += (it_third->second).second.size();
	    tag_record << '1' << endl;
	  }
	  else
	  {
	    mis[pos] += (it_third->second).first;
	    dig_mis[pos] += (it_third->second).second.size();
	    all[pos] += (it_third->second).first;
	    digital_all[pos] += (it_third->second).second.size();
	    sum_mis += (it_third->second).first;
	    sum_dig_mis += (it_third->second).second.size();
	    sum_all += (it_third->second).first;
	    sum_dig_all += (it_third->second).second.size();
	    tag_record << '0' << endl;
	  }	  
	}
      }
    }
    for (int i = 0; i != app_tag_len + 1; i++)
    {
      unformated_result << it->first << '\t' << dna_tag_map[it->first]["SUM"] << '\t';
      unformated_result << i << '\t';
      unformated_result << match[i] << '\t' << digital_match[i] << '\t';
      unformated_result << all[i] << '\t' << digital_all [i] << '\t';
      unformated_result << ((sum_match == 0) ? (0):(match[i]*100/sum_match)) << '\t' << ((sum_dig_match == 0) ? (0) : (digital_match[i]*100/sum_dig_match)) << '\t';
      unformated_result << ((sum_all == 0) ? (0):(all[i]*100/sum_all)) << '\t' << ((sum_dig_all == 0) ? (0) : (digital_all[i]*100/sum_dig_all)) << '\t';
      unformated_result << ((sum_mis == 0) ? (0) : (mis[i]*100/sum_mis)) << '\t' << ((sum_dig_mis == 0) ? (0) : (dig_mis[i]*100/sum_dig_mis)) << '\n';
    }
    unformated_result << it->first << '\t' << dna_tag_map[it->first]["SUM"] << '\t';
    unformated_result << sum_match << '\t' << sum_dig_match << '\t' << sum_mis << '\t' << sum_dig_mis << '\t' << sum_all << '\t' << sum_dig_all << '\n';
  }
  unformated_result.close();
  tag_record.close();

  run_log << "Total RNA reads: " << total_reads << endl;
  run_log << "Recovered RNA reads: " << recovered_reads << endl;
  run_log << "Recovered RNA tag diversity: " << raw_rna_tag_set.size() << endl;
  run_log << "Recovered RNA digital tag diversity: " << raw_rna_digital_tag_set.size() << endl;
  run_log << "Recovered RNA key diversity (count same tag with different lengths): " << raw_rna_key_set.size() << endl << endl;;
  run_log << "Discarded RNA reads: " << discard_reads << endl;
  run_log << "Discarded RNA reads for quality score failure: " << q_fail << endl;
  run_log << "Discarded RNA reads for template match failure: " << c_fail << endl << endl;
  run_log << "Parsed reads (counting all reads passed QC): " << parsed_reads << endl;
  double percent = (double)parsed_reads/total_reads*100;
  run_log << percent << "% of reads passed quality check\n";
  run_log << "Diversity of RNA keys passed QC: " << rna_map.size() << endl;
  run_log << "Diversity of RNA keys passed QC without matching DNA keys: " << no_dna_key_match_set.size() << endl;
  run_log << "Count of RNA keys passed QC without matching DNA keys: " << no_dna_key_match << endl << endl;
  run_log << "For all RNA reads passed QC with matching DNA keys: " << endl;
  run_log << "\t|RNA tag diversity: " << rna_tag_map.size() << endl;
  run_log << "\t|RNA key diversity: " << rna_good_key_set.size() << endl;
  run_log << "\t|RNA digital tag diversity: " << rna_digital_tag_key_div_map.size()<<endl;
  run_log << "\t|RNA reads count: " << good_reads_count << endl;
  run_log << "\t|RNA tag matched reads count: " << all_matched_reads << endl;
  run_log << "\t|RNA tag matched reads digital count: " << all_matched_digital_reads << endl << "\t|" << endl;

  vector<int> key_div_per_tag, key_count_per_tag, key_div_per_dig_tag, key_count_per_dig_tag;
  for (map<string, set<string> >::iterator it = rna_tag_map.begin(); it != rna_tag_map.end(); it++)
  {
    key_div_per_tag.push_back((it->second).size());
    key_count_per_tag.push_back(rna_key_count_map[it->first]);
  }

  for (map<string, int>::iterator it = rna_digital_tag_key_count_map.begin(); it != rna_digital_tag_key_count_map.end(); it++)
  {
    key_div_per_dig_tag.push_back(rna_digital_tag_key_div_map[it->first].size());
    key_count_per_dig_tag.push_back(rna_digital_tag_key_count_map[it->first]);
  }
  
  sort(key_div_per_tag.begin(), key_div_per_tag.end());
  sort(key_count_per_tag.begin(), key_count_per_tag.end());
  sort(key_count_per_dig_tag.begin(), key_count_per_dig_tag.end());
  sort(key_div_per_dig_tag.begin(), key_div_per_dig_tag.end());
  
  double key_div_per_tag_sum = 0, key_count_per_tag_sum = 0, key_count_per_dig_tag_sum = 0, key_div_per_dig_tag_sum = 0, key_div_per_tag_mean = 0, key_count_per_tag_mean = 0, key_div_per_dig_tag_mean = 0, key_count_per_dig_tag_mean = 0;
  double key_div_per_tag_accum = 0, key_div_per_dig_tag_accum = 0, key_count_per_tag_accum = 0, key_count_per_dig_tag_accum = 0;
  double key_div_per_tag_std = 0, key_div_per_dig_tag_std = 0, key_count_per_tag_std = 0, key_count_per_dig_tag_std = 0;

  key_div_per_tag_sum = accumulate(key_div_per_tag.begin(), key_div_per_tag.end(), 0);
  key_count_per_tag_sum = accumulate(key_count_per_tag.begin(), key_count_per_tag.end(), 0);
  key_count_per_dig_tag_sum = accumulate(key_count_per_dig_tag.begin(), key_count_per_dig_tag.end(), 0);
  key_div_per_dig_tag_sum = accumulate(key_div_per_dig_tag.begin(), key_div_per_dig_tag.end(), 0);
  
  key_div_per_tag_mean = (key_div_per_tag.size() == 0) ? (0) : (key_div_per_tag_sum / key_div_per_tag.size());
  key_count_per_tag_mean = (key_count_per_tag.size() == 0) ? (0) : (key_count_per_tag_sum/key_count_per_tag.size());
  key_div_per_dig_tag_mean = (key_div_per_dig_tag.size() == 0) ? (0) : (key_div_per_dig_tag_sum / key_div_per_dig_tag.size());
  key_count_per_dig_tag_mean = (key_count_per_dig_tag.size() == 0) ? (0) : (key_count_per_dig_tag_sum / key_count_per_dig_tag.size());
  
  for (vector<int>::iterator it = key_div_per_tag.begin(); it != key_div_per_tag.end(); it++)
  {
    key_div_per_tag_accum += (*it - key_div_per_tag_mean) * (*it - key_div_per_tag_mean);
  }
  key_div_per_tag_std = (key_div_per_tag.size() < 2) ? (0) : ( sqrt(key_div_per_tag_accum/ (key_div_per_tag.size() - 1) ) );

  for (vector<int>::iterator it = key_count_per_dig_tag.begin(); it != key_count_per_dig_tag.end(); it++)
  {
    key_count_per_dig_tag_accum += (*it - key_count_per_dig_tag_mean) * (*it - key_count_per_dig_tag_mean);
  }
  key_count_per_dig_tag_std = (key_count_per_dig_tag.size() < 2) ? (0) : ( sqrt(key_count_per_dig_tag_accum/ (key_count_per_dig_tag.size() - 1) ) );

  for (vector<int>::iterator it = key_div_per_dig_tag.begin(); it != key_div_per_dig_tag.end(); it++)
  {
    key_div_per_dig_tag_accum += (*it - key_div_per_dig_tag_mean) * (*it - key_div_per_dig_tag_mean);
  }
  key_div_per_dig_tag_std = (key_div_per_dig_tag.size() < 2) ? (0) : ( sqrt(key_div_per_dig_tag_accum/ (key_div_per_dig_tag.size() - 1) ) );

  for (vector<int>::iterator it = key_count_per_tag.begin(); it != key_count_per_tag.end(); it++)
  {
    key_count_per_tag_accum += (*it - key_count_per_tag_mean) * (*it - key_count_per_tag_mean);
  }
  key_count_per_tag_std = (key_count_per_tag.size() < 2) ? (0) : ( sqrt(key_count_per_tag_accum/ (key_count_per_tag.size() - 1) ) );

  run_log << "\t|Average key diversity per tag: " << key_div_per_tag_mean << endl;
  run_log << "\t|Median key diversity per tag: " << ( (key_div_per_tag.size() == 0) ? (0) : (key_div_per_tag[(key_div_per_tag.size()/2)]) ) << endl;
  run_log << "\t|Range key diversity per tag: " << ( (key_div_per_tag.size() == 0) ? (0) : (key_div_per_tag[key_div_per_tag.size()-1] - key_div_per_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of key diversity per tag: " << key_div_per_tag_std << endl;

  run_log << "\t|Average key count per tag: " << key_count_per_tag_mean << endl;
  run_log << "\t|Median key count per tag: " << ( (key_count_per_tag.size() == 0) ? (0) : (key_count_per_tag[(key_count_per_tag.size()/2)]) ) << endl;
  run_log << "\t|Range key count per tag: " << ( (key_count_per_tag.size() == 0) ? (0) : (key_count_per_tag[key_count_per_tag.size()-1] - key_count_per_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of key count per tag: " << key_count_per_tag_std << endl << "\t|" << endl;

  run_log << "\t|Average key diversity per digital tag: " << key_div_per_dig_tag_mean << endl;
  run_log << "\t|Median key diversity per digital tag: " << ( (key_div_per_dig_tag.size() == 0) ? (0) : (key_div_per_dig_tag[(key_div_per_dig_tag.size()/2)]) ) << endl;
  run_log << "\t|Range key diversity per digital tag: " << ( (key_div_per_dig_tag.size() == 0) ? (0) : (key_div_per_dig_tag[key_div_per_dig_tag.size()-1] - key_div_per_dig_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of key diversity per digital tag: " << key_div_per_dig_tag_std << endl;

  run_log << "\t|Average key count per digital tag: " << key_count_per_dig_tag_mean << endl;
  run_log << "\t|Median key count per digital tag: " << ( (key_count_per_dig_tag.size() == 0) ? (0) : (key_count_per_dig_tag[(key_count_per_dig_tag.size()/2)]) ) << endl;
  run_log << "\t|Range key count per digital tag: " << ( (key_count_per_dig_tag.size() == 0) ? (0) : (key_count_per_dig_tag[key_count_per_dig_tag.size()-1] - key_count_per_dig_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of key count per digital tag: " << key_count_per_dig_tag_std << endl;

  run_log.close();
  return 0;
}
