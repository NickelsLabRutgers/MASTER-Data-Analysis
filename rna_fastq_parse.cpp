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
    cerr << "Usage: " << argv[0] << " [quality score cutoff] [digital tag length] [tag length] [extra position] [key length] <output tag record file name> <run stats file name> <dna parsed file> <rna FASTQ files>\n";
    return 0;
  }

  const int qual_min = atoi(argv[1]);
  const int digital_tag_len = atoi(argv[2]), tag_len = atoi(argv[3]), extra_pos = atoi(argv[4]), key_len = atoi(argv[5]);
  const int app_tag_len = tag_len + extra_pos;
  const int total_pos = tag_len + extra_pos + 1;

  char qual_min_char;
  if ( qual_min <= -1 || qual_min >= 93 )
  {
    cerr << "Error converting qual score.\n";
    return 0;
  }
  else
  {
    qual_min_char = qual_min + 33;
  }
  
  string rna_file_name = argv[6];
  
  string sample_id = rna_file_name.substr(0,5);
  ifstream fastq_file, dna_parsed_file(argv[8]);
  ofstream tag_record(argv[6]), run_log(argv[7]);

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

  map<string, int> raw_rna_key_map, raw_rna_tag_map, raw_rna_digital_tag_map;
  map<string, pair<string, int> > dna_key_map;
  map<string, map<string, int> > dna_tag_map;
  map<string, map<string, pair<int, set<string> > > > rna_map;//{key:{tag:(count, dig_tag_set),...},...}
  map<string, set<string> > rna_tag_map;
  map<string, int> rna_key_count_map, rna_digital_tag_key_count_map;
  map<string, set<string> > rna_digital_tag_key_div_map;

  int no_dna_key_match = 0, good_reads_count = 0;
  set<string> no_dna_key_match_set, rna_good_key_set;


  while (dna_parsed_file)
  {
    dna_parsed_file >> dna_key >> dna_tag >> dna_key_count;
    if (dna_parsed_file.eof())
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
  dna_parsed_file.close();
  
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

      if(mid_index < 0 || mid_index + MID_SEQ.length() + END_SEQ.length() + key_len > bases.length() || mid_index - digital_tag_len + extra_pos < 1)
      {
        c_fail++; //Size fails.
        flag = 1;
      }
      else
      {
        rna_key = bases.substr(mid_index + MID_SEQ.length(), key_len);
        rna_dig_tag = bases.substr(0, digital_tag_len);
        rna_tag = bases.substr(digital_tag_len, mid_index - digital_tag_len + extra_pos);
        recovered_reads++;
        if (END_SEQ.compare(bases.substr(mid_index + key_len + MID_SEQ.length(), END_SEQ.length())) != 0)
        {
          flag = 1;
          c_fail++;
        }
      }

      if (flag != 1)
      {
        for(int j = 0; j < (mid_index + key_len + MID_SEQ.length() + END_SEQ.length()); j++)
        {
          if((qual.at(j) < qual_min_char) || bases.at(j) == 'N')
          {
            q_fail++; //Qual fails.
            flag = 1;
            break;
          }
        }
      }

      if(flag != 1)
      {
        rna_map[rna_key][rna_tag].first++;
        rna_map[rna_key][rna_tag].second.insert(rna_dig_tag);
        if (dna_key_map.find(rna_key) != dna_key_map.end())
        {
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
    for (set<string>::iterator itsec = (it->second).begin(); itsec != (it->second).end(); itsec++)
    {
      for (map<string, pair<int, set<string> > >::iterator it_third = (rna_map[*itsec]).begin(); it_third != (rna_map[*itsec]).end(); it_third++)
      {
        int ic_length = (it_third->first).length();
        tag_record << (it_third->first) << '\t' << (it->first) << '\t' << (tag_len + extra_pos - ic_length + 1) << '\t' << (int)((it_third->second).first) << '\t' << (int)((it_third->second).second.size()) << '\t';
        if (tag_match(it->first, it_third->first))
        {
          all_matched_reads += (it_third->second).first;
          all_matched_digital_reads += (it_third->second).second.size();
          tag_record << '1' << endl;
        }
        else
        {
          tag_record << '0' << endl;
        }
      }
    }
  }
  tag_record.close();

  run_log << "Total RNA reads: " << total_reads << endl;
  run_log << "Count of reads with unexpected structure: " << c_fail << endl;
  run_log << "Count of reads with quality score lower than cutoff: " << q_fail << endl << endl;

  run_log << "For all RNA reads passed quality checks (sequence structure and quality score): " << endl;
  run_log << "--Count:" << parsed_reads << endl;
  double percent = (double)parsed_reads/total_reads*100;
  run_log << "--Include " << percent << "% of total reads\n";
  run_log << "--Barcode diversity: " << rna_map.size() << endl;
  run_log << "--Diversity of RNA barcode without matching DNA barcode: " << no_dna_key_match_set.size() << endl;
  run_log << "--Count of RNA reads without matching DNA barcode: " << no_dna_key_match << endl << endl;
  
  run_log << "For all RNA reads passed quality checks with matching DNA barcode: " << endl;
  run_log << "\t|Count: " << good_reads_count << endl;
  run_log << "\t|TSS-region sequence diversity: " << rna_tag_map.size() << endl;
  run_log << "\t|Barcode diversity: " << rna_good_key_set.size() << endl;
  run_log << "\t|Digital tag diversity: " << rna_digital_tag_key_div_map.size()<<endl;
  run_log << "\t|Count of RNA reads with transcribed TSS-region sequence matched to DNA template: " << all_matched_reads << endl;
  run_log << "\t|Digital tag count of RNA reads with transcribed TSS-region sequence matched to DNA template: " << all_matched_digital_reads << endl << "\t|" << endl;

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

  run_log << "\t|Average barcode diversity per TSS-region: " << key_div_per_tag_mean << endl;
  run_log << "\t|Median barcode diversity per TSS-region: " << ( (key_div_per_tag.size() == 0) ? (0) : (key_div_per_tag[(key_div_per_tag.size()/2)]) ) << endl;
  run_log << "\t|Range barcode diversity per TSS-region: " << ( (key_div_per_tag.size() == 0) ? (0) : (key_div_per_tag[key_div_per_tag.size()-1] - key_div_per_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of barcode diversity per TSS-region: " << key_div_per_tag_std << endl;

  run_log << "\t|Average read count per TSS-region: " << key_count_per_tag_mean << endl;
  run_log << "\t|Median read count per TSS-region: " << ( (key_count_per_tag.size() == 0) ? (0) : (key_count_per_tag[(key_count_per_tag.size()/2)]) ) << endl;
  run_log << "\t|Range read count per TSS-region: " << ( (key_count_per_tag.size() == 0) ? (0) : (key_count_per_tag[key_count_per_tag.size()-1] - key_count_per_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of read count per TSS-region: " << key_count_per_tag_std << endl << "\t|" << endl;

  run_log << "\t|Average barcode diversity per digital tag: " << key_div_per_dig_tag_mean << endl;
  run_log << "\t|Median barcode diversity per digital tag: " << ( (key_div_per_dig_tag.size() == 0) ? (0) : (key_div_per_dig_tag[(key_div_per_dig_tag.size()/2)]) ) << endl;
  run_log << "\t|Range barcode diversity per digital tag: " << ( (key_div_per_dig_tag.size() == 0) ? (0) : (key_div_per_dig_tag[key_div_per_dig_tag.size()-1] - key_div_per_dig_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of barcode diversity per digital tag: " << key_div_per_dig_tag_std << endl;

  run_log << "\t|Average read count per digital tag: " << key_count_per_dig_tag_mean << endl;
  run_log << "\t|Median read count per digital tag: " << ( (key_count_per_dig_tag.size() == 0) ? (0) : (key_count_per_dig_tag[(key_count_per_dig_tag.size()/2)]) ) << endl;
  run_log << "\t|Range read count per digital tag: " << ( (key_count_per_dig_tag.size() == 0) ? (0) : (key_count_per_dig_tag[key_count_per_dig_tag.size()-1] - key_count_per_dig_tag[0]) ) << endl;
  run_log << "\t|Standard deviation of read count per digital tag: " << key_count_per_dig_tag_std << endl;

  run_log.close();
  return 0;
}
