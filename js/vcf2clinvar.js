
(function() 
{
  "use strict";

  var VCF2ClinVar = function() { };
  VCF2ClinVar.prototype.clinvar_report = clinvar_report;
  VCF2ClinVar.prototype.clinvar_report_json = clinvar_report_json;
  var vcf2clinvar = new VCF2ClinVar();

  if ((typeof module !== 'undefined') && (module.exports)) {
    module.exports = vcf2clinvar;
  }

  else {
    if (typeof (document) !== "undefined") window.vcf2clinvar = vcf2clinvar;
    else self['vcf2clinvar'] = vcf2clinvar;
  }


  function utf8_to_html(str) {
    var encodedStr = str.replace(/[\u00A0-\u9999<>\&]/gim, function(i) {
       return '&#'+i.charCodeAt(0)+';';
    });
    return encodedStr;
  }

  var CHROM_MAP  = {
    "1":"chr1", "chr1":"chr1", "CHR1":"chr1",
    "2":"chr2", "chr2":"chr2", "CHR2":"chr2",
    "3":"chr3", "chr3":"chr3", "CHR3":"chr3",
    "4":"chr4", "chr4":"chr4", "CHR4":"chr4",
    "5":"chr5", "chr5":"chr5", "CHR5":"chr5",
    "6":"chr6", "chr6":"chr6", "CHR6":"chr6",
    "7":"chr7", "chr7":"chr7", "CHR7":"chr7",
    "8":"chr8", "chr8":"chr8", "CHR8":"chr8",
    "9":"chr9", "chr9":"chr9", "CHR9":"chr9",
    "10":"chr10", "chr10":"chr10", "CHR10":"chr10",
    "11":"chr11", "chr11":"chr11", "CHR11":"chr11",
    "12":"chr12", "chr12":"chr12", "CHR12":"chr12",
    "13":"chr13", "chr13":"chr13", "CHR13":"chr13",
    "14":"chr14", "chr14":"chr14", "CHR14":"chr14",
    "15":"chr15", "chr15":"chr15", "CHR15":"chr15",
    "16":"chr16", "chr16":"chr16", "CHR16":"chr16",
    "17":"chr17", "chr17":"chr17", "CHR17":"chr17",
    "18":"chr18", "chr18":"chr18", "CHR18":"chr18",
    "19":"chr19", "chr19":"chr19", "CHR19":"chr19",
    "20":"chr20", "chr20":"chr20", "CHR20":"chr20",
    "21":"chr21", "chr21":"chr21", "CHR21":"chr21",
    "22":"chr22", "chr22":"chr22", "CHR22":"chr22",
    "23":"chrX", "chr23":"chrX",
    "24":"chrY", "chr24":"chrY",
    "25":"chrM", "chr25":"chrM",
    "26":"chr?", "chr26":"chr?",
    "X":"chrX", "x":"chrX", "chrx":"chrX", "CHRX":"chrX",
    "Y":"chrY", "y":"chrY", "chry":"chrY", "CHRY":"chrY",
    "MT":"chrM", "mt":"chrM", "chrM":"chrM", "CHRM":"chrM" }

  function _chrom(c) {
    if (c in CHROM_MAP) { return CHROM_MAP[c]; }
    return undefined ;
  }

  var CHROM_INDEX = {
      '1': 1, '2': 2, '3': 3, '4': 4, '5': 5,
      '6': 6, '7': 7, '8': 8, '9': 9, '10': 10,
      '11': 11, '12': 12, '13': 13, '14': 14, '15': 15,
      '16': 16, '17': 17, '18': 18, '19': 19, '20': 20,
      '21': 21, '22': 22, 'X': 23, 'Y': 24, 'M': 25, 'MT': 25,
      'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5,
      'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10,
      'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
      'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20,
      'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24, 'chrM': 25, 'chrMT': 25,
  };

  var REV_CHROM_INDEX = {
      1: 'chr1', 2: 'chr2', 3: 'chr3', 4: 'chr4', 5: 'chr5',
      6: 'chr6', 7: 'chr7', 8: 'chr8', 9: 'chr9', 10: 'chr10',
      11: 'chr11', 12: 'chr12', 13: 'chr13', 14: 'chr14', 15: 'chr15',
      16: 'chr16', 17: 'chr17', 18: 'chr18', 19: 'chr19', 20: 'chr20',
      21: 'chr21', 22: 'chr22', 23: 'chrX', 24: 'chrY', 25: 'chrMT',
  };

  vcf2clinvar.CHROM_MAP = CHROM_MAP;
  vcf2clinvar.CHROM_INDEX = CHROM_INDEX;
  vcf2clinvar.REV_CHROM_INDEX = REV_CHROM_INDEX;

  function _chrom_idx(c) {
    if (c in CHROM_INDEX) { return CHROM_INDEX[c]; }
    return 0;
  }

  function _rev_chrom_idx(c) {

    if (c in REV_CHROM_INDEX) { return REV_CHROM_INDEX[c]; }
    return "unk";

  }

  function _info_map(info_str) {
    var m = {};
    var a = info_str.split(";");
    for (var ii=0; ii<a.length; ii++) {
      var kv = a[ii].split("=");
      m[kv[0]] = kv[1];
    }
    return m;
  }

  function _parse_genotype(format_str, genotype_str) {
    var fmt_a = format_str.split(":");
    var gt_a = genotype_str.split(":");

    var idx = -1;
    for (var ii=0; ii<fmt_a.length; ii++) {
      if (fmt_a[ii] == "GT") { idx = ii; break; }
    }
    if (idx<0) { return []; }
    if (idx >= gt_a.length) { return []; }

    var gt = gt_a[idx].split(/[\|\/]/);
    for (var ii=0; ii<gt.length; ii++) {
      gt[ii] = parseInt(gt[ii]);
    }
    return gt;

  }

  function _allele_seq(ref_str, alt_str) {
    var seq_a = [ref_str];
    if (alt_str == ".") { return seq_a; }

    var alt = alt_str.split(",");
    for (var ii=0; ii<alt.length; ii++) {
      seq_a.push(alt[ii]);
    }
    return seq_a;
  }

  var CHROM_VCF_IDX = 0,
      POS_VCF_IDX = 1,
      ID_VCF_IDX = 2,
      REF_VCF_IDX = 3,
      ALT_VCF_IDX = 4,
      QUAL_VCF_IDX = 5,
      FILTER_VCF_IDX = 6,
      INFO_VCF_IDX = 7,
      FORMAT_VCF_IDX = 8,
      GT_VCF_IDX = 9;


  function _parse_frequencies(caf_str) {
    var given_freqs = caf_str.replace(/^\[/, '').replace(/\]$/, '').split(",");
    var parsed_freqs = [];

    for (var ii=0; ii<given_freqs.length; ii++) {
      if (given_freqs[ii] == ".") { parsed_freqs.push("Unknown"); }
      else { parsed_freqs.push(given_freqs[ii]); }
    }
    return parsed_freqs;
  }

  function _parse_clinvar_allele(ref_seq, alt_seqs, allele_idx, cln_idx, cln_data, frequency) {
    frequency = ((typeof frequency === "undefined") ? "Unknown" : frequency );
    var seq = ( (allele_idx==0) ? ref_seq : alt_seqs[allele_idx-1] );

    var clnsrcs = cln_data["CLNSRC"][cln_idx];
    var clnsrcids = cln_data["CLNSRCID"][cln_idx];
    var clnhgvs = cln_data["CLNHGVS"][cln_idx][0];

    var lens = [
      cln_data['CLNACC'][cln_idx].length,
      cln_data['CLNDSDB'][cln_idx].length,
      cln_data['CLNDSDBID'][cln_idx].length,
      cln_data['CLNDBN'][cln_idx].length,
      cln_data['CLNSIG'][cln_idx].length
    ];

    var n = cln_data['CLNACC'][cln_idx].length;
    for (var ii=1; ii<lens.length; ii++) {
      if (lens[ii] != n) { n=0; break; }
    }

    var records = [];
    for (var ii=0; ii<n; ii++) {

      var clndsdbs = cln_data['CLNDSDB'][cln_idx][ii].split(":");
      var clndsdbids = cln_data['CLNDSDBID'][cln_idx][ii].split(":");
      var dsdb = [];
      for (var jj=0; jj<clndsdbs.length; jj++) {
        dsdb.push( [ clndsdbs[jj], clndsdbids[jj] ] );
      }
      records.push({
        "dsdb": dsdb,
        "acc": cln_data['CLNACC'][cln_idx][ii],
        "dbn": cln_data['CLNDBN'][cln_idx][ii],
        "sig": cln_data['CLNSIG'][cln_idx][ii] });
    }

    return {
      "sequence":seq,
      "clnhgvs":clnhgvs,
      "clnsrcs":clnsrcs,
      "clnsrcids":clnsrcids,
      "records":records,
      "frequency":frequency }
  }

  function _parse_allele(ref_seq, alt_seqs, allele_idx, frequency) {
    frequency = ((typeof frequency === "undefined") ? "Unknown" : frequency );
    var seq = ( (allele_idx==0) ? ref_seq : alt_seqs[allele_idx-1] );
    return { "sequence": seq, "frequency": frequency };
  }

  function _parse_allele_data(vcf_line) {
    var vcf_fields = vcf_line.split("\t");
    var ref_seq = vcf_fields[REF_VCF_IDX];
    var alt_alleles = vcf_fields[ALT_VCF_IDX].split(",");

    var info_map = _info_map(vcf_fields[INFO_VCF_IDX]);

    var frequencies = ( ("CAF" in info_map) ? _parse_frequencies(info_map["CAF"]) : [] );

    var clnallele_keys = info_map["CLNALLE"].split(",").map(function(x){ return parseInt(x); });
    var info_clinvar_tags = ['CLNDSDB', 'CLNDSDBID', 'CLNACC', 'CLNDBN',
                             'CLNSIG', 'CLNHGVS', 'CLNSRC', 'CLNSRCID'];

    var cln_data = {};
    for (var ii in info_clinvar_tags) {
      var x = [];
      var a = info_map[info_clinvar_tags[ii]].split(",");
      for (var jj in a) { x.push( a[jj].split("|") ); }
      cln_data[info_clinvar_tags[ii]] = x;
    }
    cln_data["allele_keys"] = clnallele_keys;

    var alleles = [];
    for (var ii=0; ii<(alt_alleles.length+1); ii++) {
      var allele = {};
      var cln_idx = clnallele_keys.indexOf(ii);
      if ( (cln_idx>=0) && (frequencies.length>0) ) {
        allele = _parse_clinvar_allele(ref_seq, alt_alleles, ii, cln_idx, cln_data, frequencies[ii]);
      }

      else if (cln_idx>=0) {
        allele = _parse_clinvar_allele(ref_seq, alt_alleles, ii, cln_idx, cln_data);
      }

      else if (frequencies.length>0) {
        allele = _parse_allele(ref_seq, alt_alleles, ii, frequencies[ii]);
      }

      else {
        allele = _parse_allele(ref_seq, alt_alleles, ii);
      }

      alleles.push(allele);
    }

    return alleles;
  }

  function clinvar_report_json(clinvar, genome) {
    var cv_idx = 0, gn_idx = 0;

    var results = [];

    while ((cv_idx < clinvar.length) &&
           (gn_idx < genome.length)) {
      if (clinvar[cv_idx].length==0) { cv_idx++; continue; }
      if (genome[gn_idx].length ==0) { gn_idx++; continue; }

      if (clinvar[cv_idx][0]=='#') { cv_idx++; continue; }
      if (genome[gn_idx][0] =='#') { gn_idx++; continue; }

      var clinvar_fields = clinvar[cv_idx].split("\t");
      var genome_fields = genome[gn_idx].split("\t");

      // skip if positions don't line up
      //
      if (_chrom_idx(clinvar_fields[CHROM_VCF_IDX]) < _chrom_idx(genome_fields[CHROM_VCF_IDX])) {
      //if (CHROM_INDEX[clinvar_fields[CHROM_VCF_IDX]] < CHROM_INDEX[genome_fields[CHROM_VCF_IDX]]) {
        cv_idx++;
        continue;
      }

      else if (_chrom_idx(clinvar_fields[CHROM_VCF_IDX]) > _chrom_idx(genome_fields[CHROM_VCF_IDX])) {
      //else if (CHROM_INDEX[clinvar_fields[CHROM_VCF_IDX]] > CHROM_INDEX[genome_fields[CHROM_VCF_IDX]]) {
        gn_idx++;
        continue;
      }

      else if (parseInt(clinvar_fields[POS_VCF_IDX]) < parseInt(genome_fields[POS_VCF_IDX])) {
        cv_idx++;
        continue;
      }

      else if (parseInt(clinvar_fields[POS_VCF_IDX]) > parseInt(genome_fields[POS_VCF_IDX])) {
        gn_idx++;
        continue;
      }

      var gt = _parse_genotype(genome_fields[FORMAT_VCF_IDX], genome_fields[GT_VCF_IDX]);
      if (gt.length==0) {
        gn_idx++;
        continue;
      }

      // Only match if ref alleles match
      //
      if (clinvar_fields[REF_VCF_IDX] != genome_fields[REF_VCF_IDX]) {
        cv_idx++;
        gn_idx++;
        continue;
      }

      var gen_vcf_seq = _allele_seq(genome_fields[REF_VCF_IDX], genome_fields[ALT_VCF_IDX]);
      var genome_seq = [];
      for (var ii=0; ii<gt.length; ii++) {
        genome_seq.push( gen_vcf_seq[gt[ii]] );
      }

      var clinvar_seq = _allele_seq(clinvar_fields[REF_VCF_IDX], clinvar_fields[ALT_VCF_IDX]);

      var zygosity = "unk";
      if (gt.length==1) { zygosity = "Hem"; }
      else if (gt.length==2) {
        zygosity = ( (genome_seq[0] == genome_seq[1]) ? "Hom" : "Het" );
      }

      // Now positions have lined up, we see if there is an actual
      // sequence match.
      //

      var clinvar_allele = _parse_allele_data(clinvar[cv_idx]);

      for (var ii=0; ii<genome_seq.length; ii++) {
        if ((ii>0) && (genome_seq[ii] == genome_seq[ii-1])) { continue; }

        for (var jj=0; jj<clinvar_allele.length; jj++) {
          if (!("records" in clinvar_allele[jj])) { continue; }

          if (genome_seq[ii] == clinvar_allele[jj].sequence) {

            var name_a = [];
            for (var ri=0; ri<clinvar_allele[jj].records.length; ri++) {
              name_a.push( clinvar_allele[jj].records[ri].dbn );
            }
            var name = name_a.join(":");

            for (var ri=0; ri<clinvar_allele[jj].records.length; ri++) {

              var name = clinvar_allele[jj].records[ri].dbn;

              var match_rec = [
                genome_fields[CHROM_VCF_IDX],
                genome_fields[POS_VCF_IDX],
                clinvar_fields[REF_VCF_IDX],
                genome_seq[ii],
                name,
                clinvar_allele[jj].records[ri].sig,
                clinvar_allele[jj].frequency,
                zygosity,
                "http://www.ncbi.nlm.nih.gov/clinvar/" + clinvar_allele[jj].records[ri].acc ];

              //console.log(match_rec.join(","));
              results.push(match_rec);
            }

          }
        }
      }

      cv_idx++;
      gn_idx++;

    }

    var data = {
      "header" : [ "chrom", "pos", "ref", "alt", "name", "sig", "freq", "zyg", "url" ],
      "results": results
    };

    return data;
  }

  function clinvar_report(clinvar, genome) {
    var cv_idx = 0, gn_idx = 0;

    var results = [];

    while ((cv_idx < clinvar.length) &&
           (gn_idx < genome.length)) {
      if (clinvar[cv_idx].length==0) { cv_idx++; continue; }
      if (genome[gn_idx].length ==0) { gn_idx++; continue; }

      if (clinvar[cv_idx][0]=='#') { cv_idx++; continue; }
      if (genome[gn_idx][0] =='#') { gn_idx++; continue; }

      var clinvar_fields = clinvar[cv_idx].split("\t");
      var genome_fields = genome[gn_idx].split("\t");

      // skip if positions don't line up
      //
      if (_chrom_idx(clinvar_fields[CHROM_VCF_IDX]) < _chrom_idx(genome_fields[CHROM_VCF_IDX])) {
      //if (CHROM_INDEX[clinvar_fields[CHROM_VCF_IDX]] < CHROM_INDEX[genome_fields[CHROM_VCF_IDX]]) {
        cv_idx++;
        continue;
      }

      else if (_chrom_idx(clinvar_fields[CHROM_VCF_IDX]) > _chrom_idx(genome_fields[CHROM_VCF_IDX])) {
      //else if (CHROM_INDEX[clinvar_fields[CHROM_VCF_IDX]] > CHROM_INDEX[genome_fields[CHROM_VCF_IDX]]) {
        gn_idx++;
        continue;
      }

      else if (parseInt(clinvar_fields[POS_VCF_IDX]) < parseInt(genome_fields[POS_VCF_IDX])) {
        cv_idx++;
        continue;
      }

      else if (parseInt(clinvar_fields[POS_VCF_IDX]) > parseInt(genome_fields[POS_VCF_IDX])) {
        gn_idx++;
        continue;
      }

      var gt = _parse_genotype(genome_fields[FORMAT_VCF_IDX], genome_fields[GT_VCF_IDX]);
      if (gt.length==0) {
        gn_idx++;
        continue;
      }

      // Only match if ref alleles match
      //
      if (clinvar_fields[REF_VCF_IDX] != genome_fields[REF_VCF_IDX]) {
        cv_idx++;
        gn_idx++;
        continue;
      }

      var gen_vcf_seq = _allele_seq(genome_fields[REF_VCF_IDX], genome_fields[ALT_VCF_IDX]);
      var genome_seq = [];
      for (var ii=0; ii<gt.length; ii++) {
        genome_seq.push( gen_vcf_seq[gt[ii]] );
      }

      var clinvar_seq = _allele_seq(clinvar_fields[REF_VCF_IDX], clinvar_fields[ALT_VCF_IDX]);

      var zygosity = "unk";
      if (gt.length==1) { zygosity = "Hem"; }
      else if (gt.length==2) {
        zygosity = ( (genome_seq[0] == genome_seq[1]) ? "Hom" : "Het" );
      }

      // Now positions have lined up, we see if there is an actual
      // sequence match.
      //

      var clinvar_allele = _parse_allele_data(clinvar[cv_idx]);

      for (var ii=0; ii<genome_seq.length; ii++) {
        if ((ii>0) && (genome_seq[ii] == genome_seq[ii-1])) { continue; }

        for (var jj=0; jj<clinvar_allele.length; jj++) {
          if (!("records" in clinvar_allele[jj])) { continue; }

          if (genome_seq[ii] == clinvar_allele[jj].sequence) {

            var name_a = [];
            for (var ri=0; ri<clinvar_allele[jj].records.length; ri++) {
              name_a.push( clinvar_allele[jj].records[ri].dbn );
            }
            var name = name_a.join(":");

            for (var ri=0; ri<clinvar_allele[jj].records.length; ri++) {

              var name = clinvar_allele[jj].records[ri].dbn;

              var match_rec = [
                genome_fields[CHROM_VCF_IDX],
                genome_fields[POS_VCF_IDX],
                clinvar_fields[REF_VCF_IDX],
                genome_seq[ii],
                name,
                clinvar_allele[jj].records[ri].sig,
                clinvar_allele[jj].frequency,
                zygosity,
                "http://www.ncbi.nlm.nih.gov/clinvar/" + clinvar_allele[jj].records[ri].acc ];

              //console.log(match_rec.join(","));
              results.push(match_rec.join(","));
            }

          }
        }
      }

      cv_idx++;
      gn_idx++;

    }

    return results;
  }

})();           
