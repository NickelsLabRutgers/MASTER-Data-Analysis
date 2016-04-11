/*-------------------------------------------------------------------------------
  Input: fastq data
  Output: (1) parsed files containing sequences passed quality control:
  '7nts after BEG_SEQ'\t'15nts after MID_SEQ' endl
  (2) discarded files containing sequences failed quality control:
  'complete sequence' endl
  -------------------------------------------------------------------------------*/
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <math.h>
#include <utility>

using namespace std;

string convertInt(int number)
{
  stringstream ss;
  ss << number;
  return ss.str();
}

int sum_map_func(int x, pair<string, int> p) 
{
  //Helper function for summing up read counts in {read:count} map.
  return x + p.second;
}

bool par_map_max_comp_func(pair<string, int> i, pair<string, int> j)
{
  //Helper function for getting max read count in {read:count} map.
  return i.second < j.second;
}

/*-----------------Reads Structure-------------------------------------------------------------
 * NN...N-[BEG_SEQ]-[N7]-[MID_SEQ]-[N15]-[END_SEQ]-N...NNN
 * BEG_SEQ, MID_SEQ, END_SEQ: Fixed sequences used for matching
 * N7: tag, 7 or 10 nucleotides long random transcription start region
 * N15: key, 15 nucleotides long random key region used for mapping RNA transcripts to DNA templates.
 ---------------------------------------------------------------------------------------------*/
const string BEG_SEQ="AGGCTTGACACTTTATGCTTCGGCTCGTATAATGTG"; //36
const string MID_SEQ="GTGAGCGGATAACAAT"; //16
const string END_SEQ="TGGAA"; //5
const string QUAL_SCORE = " !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

int main(int argc, char * argv[])
{
  if(argc < 8)
  {
    /*-----------Arguments-------------------
     argv[1]: Sequencing quality score for cutoff. Integer.
     argv[2]: tag(random TSS region) length.
     argv[3]: extra position number considered downstream tag region.
     argv[4]: key length.
     argv[5]: File name for output DNA template file.
     argv[6]: File name for output stats file.
     argv[7]: Input DNA sequencing library raw data files.
     */
    cerr << "Usage: " << argv[0] << " [quality score cutoff] [tag length] [extra position] [key length] <output DNA parsed file name> <output DNA stats file name> <FASTQ files>\n";
    return 0;
  }

  const int qual_min = atoi(argv[1]), tag_len = atoi(argv[2]), extra_pos = atoi(argv[3]), key_len = atoi(argv[4]);
  char qual_min_char;
  if ( qual_min <= -1 || qual_min >= 95 )
  {
    cerr << "Error converting qual score.\n";
    return 0;
  }
  else
  {
    qual_min_char = QUAL_SCORE[qual_min];
  }
  
  string parsed_out = argv[5];
  string stats_out = argv[6];
  string file_name;
  ifstream fastq_file;
  ofstream parsed_file, run_log;
  parsed_file.open(parsed_out.c_str());
  run_log.open(stats_out.c_str());

  int flag = 0;
  int s_fail_c = 0;
  int q_fail_c = 0;
  int total_reads = 0;
  int parsed_reads = 0;
  int multi_tag_discarded_reads_diversity = 0;
  int multi_tag_discarded_reads_count = 0;
  int tag_count = 0;
  int key_count = 0;
  double percent = 0;
  double key_div_sum = 0, key_count_sum = 0, key_div_mean = 0, key_count_mean = 0, key_count_accum = 0, key_div_accum = 0, key_count_std = 0, key_div_std = 0;
  string key, tag;
  map<string, map<string, int> > parsed_map;
  map<string, map<string, int> > tag_map;
  vector<int> keys_per_tag_div, keys_per_tag_count;
  map<string, int> raw_key_count_map, raw_tag_count_map;
  set<string> raw_key_set, raw_tag_set, raw_digital_tag_set;
  string p_name = argv[5];

  int begin_pos = 0, mid_start = 0, end_start = 0;

  for(int i = 7; i < argc; i++)
  {
    file_name = argv[i];
    fastq_file.open(file_name.c_str());
    if(!fastq_file)
    {
      cerr << "Error reading file\n";
      return 0;
    }
    string label_hold,bases,pos,qual;
    while(fastq_file)
    {
      //FastQ data: (1) sequence identifier (2) sequence (3) + followed by same sequence identifier (4) quality values
      getline(fastq_file,label_hold);
      getline(fastq_file,bases);
      getline(fastq_file,pos);
      getline(fastq_file,qual);
      if(fastq_file.eof())
        break;

      total_reads++;

      flag = 0;

      //BEG_SEQ-7 nts-MID_SEQ-15 nts-END_SEQ
      //         Tag            Key
      

      begin_pos = bases.find(BEG_SEQ);
      if (begin_pos != -1 && begin_pos + BEG_SEQ.length() + MID_SEQ.length() + END_SEQ.length() + key_len + tag_len <= bases.length() )
      {
        mid_start = begin_pos + BEG_SEQ.length() + tag_len;
        end_start = mid_start + MID_SEQ.length() + key_len;
        key = bases.substr(mid_start + MID_SEQ.length(), key_len); //key
        tag = bases.substr(begin_pos + BEG_SEQ.length(), tag_len + extra_pos); //tag
      }
      else
      {
        flag = 1;
        s_fail_c++;
        begin_pos = 0;
        mid_start = BEG_SEQ.length() + tag_len;
        end_start = mid_start + MID_SEQ.length() + key_len;
      }

      if (flag != 1) 
      {
        if (MID_SEQ.substr(extra_pos, MID_SEQ.length() - extra_pos).compare(bases.substr(mid_start + extra_pos, MID_SEQ.length() - extra_pos)) != 0 || END_SEQ.compare(bases.substr(end_start, END_SEQ.length())) != 0)
        {
          s_fail_c++;
          flag = 1;
        }
      }

      if (flag != 1) 
      {
        for (int j = begin_pos; j < begin_pos + BEG_SEQ.length() + MID_SEQ.length() + END_SEQ.length() + key_len + tag_len; j++)
        {
          if((qual.at(j) < qual_min_char) || bases.at(j) == 'N')
          {
            q_fail_c++;
            flag = 1;
            break;
          }
        }
      }

      if(flag != 1)
      {
        raw_key_set.insert(key);
        raw_tag_set.insert(tag);
        if (raw_key_count_map.find(key) == raw_key_count_map.end())
        {
          raw_key_count_map[key] = 1;
        }
        else
        {
          raw_key_count_map[key]++;
        }

        if (raw_tag_count_map.find(tag) == raw_tag_count_map.end())
        {
          raw_tag_count_map[tag] = 1;
        }
        else
        {
          raw_tag_count_map[tag]++;
        }

        if (parsed_map.find(key) != parsed_map.end())
        {
          //parsed_map = {key: {tag:read_count} ... }
          if (parsed_map[key].find(tag) != parsed_map[key].end())
            parsed_map[key][tag]++;
          else
            parsed_map[key][tag] = 1;
        }
        else
        {
          parsed_map[key][tag] = 1;
        }
      }
    }
    fastq_file.close();
  }

  for (map<string, map<string, int> >::iterator it = parsed_map.begin(); it != parsed_map.end(); it++) 
  {
    //Iterate through passed map, if key mapped to different tags and percentage of max_tag_count is below 0.9, discard key.
    //If max tag read count of a key was lower than 10, discarde key. Otherwise insert into tag_map.
    //tag_map = {tag:{key:count} ... } keys contain "SUM", which is the sum of read counts of tag.
    int max_count = 0, sum_count = 0;
    max_count = max_element((it->second).begin(), (it->second).end(), par_map_max_comp_func) -> second;
    sum_count = accumulate((it->second).begin(), (it->second).end(), 0, sum_map_func);
    double key_percentage = 0;
    key_percentage = static_cast<double>(max_count) / static_cast<double>(sum_count);
    if ( key_percentage < 0.9 ) 
    {
      for (map<string, int>::iterator itn = (it->second).begin(); itn != (it->second).end(); itn++) 
      {
        multi_tag_discarded_reads_diversity++;
        multi_tag_discarded_reads_count += itn->second;
      }
    } 
    else if ( max_count >= 10 )
    {
      for (map<string, int>::iterator itn = (it->second).begin(); itn != (it->second).end(); itn++) 
      {
        if ((itn->second) == max_count)
        {
          tag_map[itn->first][it->first] = itn->second;
          if (tag_map[itn->first].find("SUM") == tag_map[itn->first].end()) {
            tag_map[itn->first]["SUM"] = itn->second;
          } 
          else 
          {
            tag_map[itn->first]["SUM"] += itn->second;
          }
        } 
        else 
        {
           multi_tag_discarded_reads_diversity++;
           multi_tag_discarded_reads_count += itn->second;
        }
      }
    }
  }
  
  for (map<string, map<string, int> >::iterator it = tag_map.begin(); it != tag_map.end(); it++) 
  {
    for (map<string, int>::iterator itn = (it->second).begin(); itn != (it->second).end(); itn++) 
    {
      if (itn->first != "SUM") 
      {
        parsed_file << (itn->first) << '\t' << (it->first) << '\t' << (itn->second) << '\n';
        parsed_reads += (itn->second);
        key_count++;
      }
    }
    tag_count++;
    keys_per_tag_div.push_back((it->second).size() - 1);
    keys_per_tag_count.push_back((it->second)["SUM"]);
  }

  sort(keys_per_tag_div.begin(), keys_per_tag_div.end());
  sort(keys_per_tag_count.begin(), keys_per_tag_count.end());
  
  run_log << "Total reads: " << total_reads << endl;
  run_log << "Count of reads with unexpected structure: " << s_fail_c << endl;
  run_log << "Count of reads with quality score lower than cutoff: " << q_fail_c << endl << endl;
  run_log << "For all DNA reads passed quality checks (sequence structure and quality score):" << endl;
  run_log << "--Count: " << (accumulate(raw_key_count_map.begin(), raw_key_count_map.end(), 0, sum_map_func)) << endl;
  run_log << "--Barcode diversity: " << raw_key_set.size() << endl;
  run_log << "--N7 diversity: " << raw_tag_set.size() << endl;
  run_log << "--Diversity of barcodes mapped to multiple N7s: " << multi_tag_discarded_reads_diversity << '\n';
  run_log << "--Count of reads with barcodes mapped to multiple N7s: " << multi_tag_discarded_reads_count << '\n' << '\n';
  
  run_log << "For all DNA reads passed quality, multiple N7s mapping, and read count checks:" << endl;
  run_log << "--Count: " << parsed_reads << endl;
  percent = static_cast<double> (parsed_reads)/total_reads*100;
  run_log << "--Include " << percent << "% of total reads.\n";
  run_log << "--N7 diversity (passed all checks): " << tag_count << endl;
  run_log << "--Barcode diversity (passed all check): " << key_count << endl;

  

  key_div_sum = accumulate(keys_per_tag_div.begin(), keys_per_tag_div.end(), 0);
  key_div_mean = (keys_per_tag_div.size() == 0) ? (0) : (key_div_sum / keys_per_tag_div.size());
  for (vector<int>::iterator it = keys_per_tag_div.begin(); it != keys_per_tag_div.end(); it++)
  {
    key_div_accum += (*it - key_div_mean) * (*it - key_div_mean);
  }
  
  key_div_std = (keys_per_tag_div.size() < 2) ? (0) : ( sqrt(key_div_accum/ (keys_per_tag_div.size() - 1) ) );

  key_count_sum = accumulate(keys_per_tag_count.begin(), keys_per_tag_count.end(), 0);
  key_count_mean = (keys_per_tag_count.size() == 0) ? (0) : (key_count_sum / keys_per_tag_count.size());
  for (vector<int>::iterator it = keys_per_tag_count.begin(); it != keys_per_tag_count.end(); it++)
  {
    key_count_accum += (*it - key_count_mean) * (*it - key_count_mean);
  }
  key_count_std =(keys_per_tag_count.size() < 2) ? (0) : ( sqrt(key_count_accum/ (keys_per_tag_count.size() - 1) ) );

  run_log << "--Average barcode diversity per N7: " << key_div_mean << endl;
  run_log << "--Median barcode diversity per N7: " << ( (keys_per_tag_div.size() == 0) ? (0) : (keys_per_tag_div[(keys_per_tag_div.size()/2)]) ) << endl;
  run_log << "--Range of barcode diversity per N7: " << ( (keys_per_tag_div.size() == 0) ? (0) : (keys_per_tag_div[keys_per_tag_div.size()-1] - keys_per_tag_div[0]) ) << endl;
  run_log << "--Standard deviation of barcode diversity per N7: " << key_div_std << endl;
  run_log << "--Average counts per N7: " << key_count_mean << endl;
  run_log << "--Median counts per N7: " << ( (keys_per_tag_count.size() == 0) ? (0) : (keys_per_tag_count[(keys_per_tag_count.size()/2)]) ) << endl;
  run_log << "--Range of counts per N7: " << ( (keys_per_tag_count.size() == 0) ? (0) : (keys_per_tag_count[keys_per_tag_count.size()-1] - keys_per_tag_count[0]) ) << endl;
  run_log << "--Standard deviation of counts per N7: " << key_count_std << endl;
  
  fastq_file.close();
  parsed_file.close();
  run_log.close();
  return 0;
}
