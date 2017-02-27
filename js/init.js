console.log("in init");

var g_clinvar_ready = false;
var g_genome_ready = false;

var g_clinvar_lines = undefined;
var g_genome_lines = undefined;

function load_file(e) {
  var data = e.target.result;
  //var data_str = data;
  //var data_str = atob(data);

  var data_parts = data.split(",");
  var data_str = atob(data_parts[1]);

  console.log("got:", data_str.length);

  g_genome_lines = data_str.split("\n");
  report_semaphore();
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

  //var data = clinvar_report_json(g_clinvar_lines, g_genome_lines);
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
}
