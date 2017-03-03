console.log("in init");

var g_clinvar_ready = false;
var g_genome_ready = false;

var g_clinvar_lines = undefined;
var g_genome_lines = undefined;

var g_23andme_ref_b37 = {};
var g_ancestrydna_ref_b37 = {};

function load_gt_reference(ref_map, ref_lines) {
  for (var i=0; i<ref_lines.length; i++) {
    var field = ref_lines[i].split("\t");
    ref_map[ field[0] + "\t" + field[1] ] = field[2];
  }
}

function convert_23andme_to_vcf(ref_map, lines) {
  var fields = [];
  var key = "";
  var vcf_lines = [];

  vcf_lines.push("##fileformat=VCFv4.2");
  vcf_lines.push('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">');
  vcf_lines.push(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GENOTYPE"].join("\t"));

  for (var i=0; i<lines.length; i++) {
    if (lines[i].length==0) { continue; }
    if (lines[i][0] == '#') { continue; }

    //fields = lines[i].split("\t");
    fields = lines[i].replace(/(\r\n|\r|\n)/gm, '').split("\t");
    key = fields[1] + "\t" + fields[2];

    // skip no calls
    //
    if ((fields[3] == "--") || 
        (fields[3] == "-")) { continue; }

    //var gt = [ fields[3][0], fields[3][1] ];
    var gt = [ fields[3][0] ];
    if (fields[3].length>1) { gt.push(fields[3][1]); }

    var gt_map = {};
    gt_map[ ref_map[key] ] = 0;
    var gt_idx = [];
    var alt_str = "";
    var count = 1;
    for (var idx=0; idx<gt.length; idx++) {
      if (gt[idx] in gt_map) {
        gt_idx.push(gt_map[gt[idx]]);
        continue;
      }
      if (alt_str.length>0) { alt_str += ","; }
      alt_str += gt[idx];
      gt_map[gt[idx]] = count;
      gt_idx.push(count);
      count += 1;
    }

    //if (i<10) { console.log("debug", fields, key, ref_map[key], gt, gt_map, gt_idx, alt_str, count); }

    var chrom = vcf2clinvar.CHROM_MAP[fields[1]];
    if (alt_str.length==0) { alt_str = "."; }

    if (key in ref_map) {
      vcf_lines.push(
          //fields[1] +                 // chrom
          chrom +                 // chrom
          "\t" + fields[2] +          // pos
          "\t" + fields[0] +          // id
          "\t" + ref_map[key] +       // ref
          "\t" + alt_str +            // alt
          "\t" + "." +                // qual
          "\t" + "." +                // filter
          "\t" + "." +                // info
          "\t" + "GT" +               // format
          "\t" + gt_idx.join("/"));   // genotype
    }
  }

  //for (var i=0; i<10; i++) { console.log(i, vcf_lines[i]); }

  return vcf_lines;
}

function convert_ancestrydna_to_vcf(ref_map, lines) {
  var fields = [];
  var key = "";
  var vcf_lines = [];

  vcf_lines.push("##fileformat=VCFv4.2");
  vcf_lines.push('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">');
  vcf_lines.push(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GENOTYPE"].join("\t"));

  for (var i=0; i<lines.length; i++) {
    if (lines[i].length==0) { continue; }
    if (lines[i][0] == '#') { continue; }

    //fields = lines[i].split("\t");
    fields = lines[i].replace(/(\r\n|\r|\n)/gm, '').split("\t");
    key = fields[1] + "\t" + fields[2];

    // skip no calls
    //
    if ((fields[3] == "--") || 
        (fields[3] == "-")  ||
        (fields[3] == ".")  ||
        (fields[4] == "--") ||
        (fields[4] == "-")  ||
        (fields[4] == ".")) { continue; }

    //var gt = [ fields[3][0], fields[3][1] ];
    var gt = [ fields[3], fields[4] ];
    //if (fields[3].length>1) { gt.push(fields[3][1]); }

    var gt_map = {};
    gt_map[ ref_map[key] ] = 0;
    var gt_idx = [];
    var alt_str = "";
    var count = 1;
    for (var idx=0; idx<gt.length; idx++) {
      if (gt[idx] in gt_map) {
        gt_idx.push(gt_map[gt[idx]]);
        continue;
      }
      if (alt_str.length>0) { alt_str += ","; }
      alt_str += gt[idx];
      gt_map[gt[idx]] = count;
      gt_idx.push(count);
      count += 1;
    }

    //if (i<10) { console.log("debug", fields, key, ref_map[key], gt, gt_map, gt_idx, alt_str, count); }

    var chrom = vcf2clinvar.CHROM_MAP[fields[1]];
    if (alt_str.length==0) { alt_str = "."; }

    if (key in ref_map) {
      vcf_lines.push(
          //fields[1] +                 // chrom
          chrom +                 // chrom
          "\t" + fields[2] +          // pos
          "\t" + fields[0] +          // id
          "\t" + ref_map[key] +       // ref
          "\t" + alt_str +            // alt
          "\t" + "." +                // qual
          "\t" + "." +                // filter
          "\t" + "." +                // info
          "\t" + "GT" +               // format
          "\t" + gt_idx.join("/"));   // genotype
    }
  }

  //for (var i=0; i<10; i++) { console.log(i, vcf_lines[i]); }

  return vcf_lines;
}

function load_file(e) {
  var data = e.target.result;
  //var data_str = data;
  //var data_str = atob(data);

  var data_parts = data.split(",");
  var data_str = atob(data_parts[1]);

  console.log("got:", data_str.length);

  //g_genome_lines = data_str.split("\n");
  var data_lines = data_str.split("\n");
  if (data_lines.length==0) { return; }
  if (data_lines[0].search(/##fileformat=VCF/)>=0) {
    g_genome_lines = data_lines;
    report_semaphore();
    return;
  }

  if (data_lines[0].search(/23and[mM]e/)>=0) {
    console.log(">>> 23andme data?");

    //for (var i=0; i<10; i++) { console.log("inp:", data_lines[i]); }

    var vcf_lines = convert_23andme_to_vcf(g_23andme_ref_b37, data_lines);

    //for (var i=0; i<5; i++) { console.log(vcf_lines[i]); }
    //for (var i=(vcf_lines.length-5); i<vcf_lines.length; i++) { console.log(vcf_lines[i]); }


    return;
  }

  if (data_lines[0].search(/Ancestry[dD][nN][aA]/)>=0) {
    console.log(">>> Ancestry DNA?");

    var vcf_lines = convert_ancestrydna_to_vcf(g_ancestrydna_ref_b37, data_lines);

    for (var i=0; i<5; i++) { console.log(vcf_lines[i]); }
    for (var i=(vcf_lines.length-5); i<vcf_lines.length; i++) { console.log(vcf_lines[i]); }


    return;
  }


  console.log("unknown format");
}

function process_genome_upload() {
  var inp_id = "genome_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(fn);

  var reader = new FileReader();
  //reader.onload = function(x) { console.log("file ul, got:", x.target.result); };
  reader.onload = load_file;
  reader.readAsDataURL(fn);
}

function remove_spinner() {
  $("#spinner").css("display", "none");
  $("#main").css("display", "inline");

}

function report_semaphore() {
  if ((!g_genome_ready) || (!g_clinvar_ready)) {
    return;
  }

  var _sig_lookup = {
    "" : "all",
    "0": "unknown",
    "1": "untested",
    "2": "non-pathogenic",
    "3": "probably non-pathogenic",
    "4": "probably pathogenic",
    "5": "pathogenic",
    "6": "affecting drug response",
    "7": "affecting histocompatibility",
    "255": "other"
  };


  /*
  var res = clinvar_report(g_clinvar_lines, g_genome_lines);
  console.log("result:", res.length, "resulting lines");
  for (var ii=0; ii<res.length; ii++) {
    console.log(res[ii]);
    if (ii>10) { break; }
  }
  */

  var data = vcf2clinvar.clinvar_report_json(g_clinvar_lines, g_genome_lines);
  var res = data.results;

  console.log(data);

  var variants = [];
  for (var ii=0; ii<res.length; ii++) {
    var v = {
			"chrom": res[ii][0],
			"pos": res[ii][1],
			"ref_allele": res[ii][2],
			"alt_allele": res[ii][3],
			"name": res[ii][4],
			"clinical_significance": res[ii][5],
			"clinical_significance_descr": _sig_lookup[res[ii][5]],
			"allele_freq": res[ii][6],
			"zygosity": res[ii][7],
			"acc_url": res[ii][8],
			"html_link_name" : "<a href='" + res[ii][8] + "'>" + res[ii][4] + "</a>"
    };
	  variants.push(v);

    //if (ii>10) { break; }
  }

  var $table = $('table');
  $(function() {
    $table.bootstrapTable({data:variants});

    // Do the default filter
    //
    $table.bootstrapTable('filterBy', {clinical_significance_descr: ["pathogenic"]});

    // kind of hacky, but seems to work...set the displayed default filter to 'pathogenic'
    //
    $(".clinical_significance_descr").val("pathogenic");

    // finally sort by allele_freq
    //
    $("th[data-field='allele_freq']").children(".sortable").click().click();

  });

  remove_spinner();
}

function on_ready() {
  var cln_wurk = new Worker("js/clinvar-worker.js");
  cln_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    g_clinvar_lines = raw_data.split("\n");

    g_clinvar_ready = true;

    report_semaphore();
    console.log("clinvar:", g_clinvar_lines.length, "lines loaded");
  });

  var example_wurk = new Worker("js/example-vcf-worker.js");
  example_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    g_genome_lines = raw_data.split("\n");

    g_genome_ready = true;

    report_semaphore();

    console.log("example:", g_genome_lines.length, "lines loaded");
  });

  var ttam_ref_wurk = new Worker("js/ref-b37-worker.js");
  ttam_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_23andme_ref_b37, raw_lines);
    console.log("g_23andme_ref loaded");
  });
  ttam_ref_wurk.postMessage("../data/23andme_reference_b37.txt.gz");

  var adna_ref_wurk = new Worker("js/ref-b37-worker.js");
  adna_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_ancestrydna_ref_b37, raw_lines);
    console.log("g_ancestrydna_ref loaded");
  });
  adna_ref_wurk.postMessage("../data/ancestrydna_reference_b37.txt.gz");

}
