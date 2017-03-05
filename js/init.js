console.log("in init");

var g_clinvar_ready = false;
var g_clinvar_queued = false;
var g_genome_ready = false;

var g_clinvar_lines = undefined;
var g_genome_lines = undefined;
var g_data_lines = undefined;

var g_ref_23andme_ready  = false;
var g_ref_23andme_queued = false;

var g_ref_ancestrydna_ready  = false;
var g_ref_ancestrydna_queued = false;

var g_ref_23andme_b37 = {};
var g_ref_ancestrydna_b37 = {};

var g_process_data = false;

function load_gt_reference(ref_map, ref_lines) {
  for (var i=0; i<ref_lines.length; i++) {
    var field = ref_lines[i].split("\t");
    ref_map[ field[0] + "\t" + field[1] ] = field[2];
  }
}

function on_ready() {
  init_table();
}

//                            _          
//  __ ___ _ ___ _____ _ _ __(_)___ _ _  
// / _/ _ \ ' \ V / -_) '_(_-< / _ \ ' \ 
// \__\___/_||_\_/\___|_| /__/_\___/_||_|
//                                       


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

    fields = lines[i].replace(/(\r\n|\r|\n)/gm, '').split("\t");
    key = fields[1] + "\t" + fields[2];

    // skip no calls
    //
    if ((fields[3] == "--") || 
        (fields[3] == "-")) { continue; }

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

// Notes: (and TODO)
// chromosomes 25 are 'PARA'
// and I think 26 is chrM
//
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

    var gt = [ fields[3], fields[4] ];

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

    var chrom = vcf2clinvar.CHROM_MAP[fields[1]];

    //FOR NOW
    //ignore 25, map 26 to chrM
    //

    //if (fields[1] == "25") { chrom = "PARA"; }
    if (fields[1] == "25") { continue; }
    if (fields[1] == "26") { chrom = "chrM"; }


    if (alt_str.length==0) { alt_str = "."; }

    if (key in ref_map) {
      vcf_lines.push(
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

  return vcf_lines;
}

//  _              _ _           
// | |___  __ _ __| (_)_ _  __ _ 
// | / _ \/ _` / _` | | ' \/ _` |
// |_\___/\__,_\__,_|_|_||_\__, |
//                         |___/ 

function auto_io(data) {
  var data_parts = data.split(",");

  console.log(">>>", data_parts[0]);

  var data_str = atob(data_parts[1]);

  if (data_parts[0] == "data:application/gzip;base64") {

    console.log("... gunzipping");

    // contortions to get it into a format for gunzip 
    //
    var uInt8Array = new Uint8Array(new ArrayBuffer(data_str.length));
    for (var ii=0; ii<data_str.length; ii++) {
      uInt8Array[ii] = data_str.charCodeAt(ii);
    }

    var gunzip = new Zlib.Gunzip(uInt8Array);
    var unpacked_uInt8Array = gunzip.decompress();
    data_str = new TextDecoder("utf-8").decode(unpacked_uInt8Array);

    var xx = data_str.slice(0,100);
    console.log(">>>", xx.length, xx);
  }

  else if (data_parts[0] == "data:application/zip;base64") {

    console.log("... unzipping");

    // contortions to get it into a format for gunzip 
    //
    var uInt8Array = new Uint8Array(new ArrayBuffer(data_str.length));
    for (var ii=0; ii<data_str.length; ii++) {
      uInt8Array[ii] = data_str.charCodeAt(ii);
    }

    // Take the first file in the zip archive
    //
    var unzip = new Zlib.Unzip(uInt8Array);
    var filenames = unzip.getFilenames();
    var unpacked_uInt8Array = unzip.decompress(filenames[0]);
    data_str = new TextDecoder("utf-8").decode(unpacked_uInt8Array);

    var xx = data_str.slice(0,100);
    console.log(">>>", xx.length, xx);

  }

  return data_str;
}

function load_clinvar_file(e) {
  var data_str = auto_io(e.target.result);
  g_clinvar_lines = data_str.split("\n");
  g_clinvar_ready = true;
  g_clinvar_queued = false;
  report_semaphore();

  if (g_clinvar_lines.length>2) {
    $("#clinvar_vcf_info").text( g_clinvar_lines[1] );
  }
  else {
    $("#clinvar_vcf_info").text( $("#clinvar_upload").val() );
  }

}

function load_ref_23andme_file(e) {
  var data_str = auto_io(e.target.result);
  var raw_lines = data_str.split("\n");
  g_ref_23andme_b37 = {};
  load_gt_reference(g_ref_23andme_b37, raw_lines);

  g_ref_23andme_ready = true;
  g_ref_23andme_queued = false;


  data_23andme_semaphore();

  $("#ref_23andme_info").text( $("#ref_23andme_upload").val() );

}

function load_ref_ancestrydna_file(e) {
  var data_str = auto_io(e.target.result);
  var raw_lines = data_str.split("\n");
  g_ref_ancestrydna_b37 = {};
  load_gt_reference(g_ref_ancestrydna_b37, raw_lines);

  g_ref_ancestrydna_ready = true;
  g_ref_ancestrydna_queued = false;

  data_ancestrydna_semaphore();

  $("#ref_ancestrydna_info").text( $("#ref_ancestrydna_upload").val() );

}

function load_genome_file(e) {

  var data_str = auto_io(e.target.result);

  //DEBUG
  console.log("load_genome_file...");
  console.log("got...", data_str.slice(0, 100));

  /*
  var data = e.target.result;
  var data_parts = data.split(",");

  console.log(">>>", data_parts[0]);

  var data_str = atob(data_parts[1]);

  if (data_parts[0] == "data:application/gzip;base64") {

    console.log("... gunzipping");

    // contortions to get it into a format for gunzip 
    //
    var uInt8Array = new Uint8Array(new ArrayBuffer(data_str.length));
    for (var ii=0; ii<data_str.length; ii++) {
      uInt8Array[ii] = data_str.charCodeAt(ii);
    }

    var gunzip = new Zlib.Gunzip(uInt8Array);
    var unpacked_uInt8Array = gunzip.decompress();
    data_str = new TextDecoder("utf-8").decode(unpacked_uInt8Array);

    var xx = data_str.slice(0,100);
    console.log(">>>", xx.length, xx);
  }

  else if (data_parts[0] == "data:application/zip;base64") {

    console.log("... unzipping");

    // contortions to get it into a format for gunzip 
    //
    var uInt8Array = new Uint8Array(new ArrayBuffer(data_str.length));
    for (var ii=0; ii<data_str.length; ii++) {
      uInt8Array[ii] = data_str.charCodeAt(ii);
    }

    // Take the first file in the zip archive
    //
    var unzip = new Zlib.Unzip(uInt8Array);
    var filenames = unzip.getFilenames();
    var unpacked_uInt8Array = unzip.decompress(filenames[0]);
    data_str = new TextDecoder("utf-8").decode(unpacked_uInt8Array);

    var xx = data_str.slice(0,100);
    console.log(">>>", xx.length, xx);


  }
  */

  g_genome_lines = undefined;

  var data_lines = data_str.split("\n");
  if (data_lines.length==0) { return; }
  if (data_lines[0].search(/##fileformat=VCF/)>=0) {

    console.log("vcf genome file...");

    g_genome_lines = data_lines;
    g_genome_ready = true;

    report_semaphore();
    return;
  }

  if (data_lines[0].search(/23and[mM]e/)>=0) {
    console.log(">>> 23andme data?");

    g_data_lines = data_lines;
    g_process_data = true;

    data_23andme_semaphore();

    //g_genome_lines = convert_23andme_to_vcf(g_ref_23andme_b37, data_lines);
    //report_semaphore();
    return;
  }

  if (data_lines[0].search(/Ancestry[dD][nN][aA]/)>=0) {
    console.log(">>> Ancestry DNA?");

    g_data_lines = data_lines;
    g_process_data = true;
    data_ancestrydna_semaphore();

    //g_genome_lines = convert_ancestrydna_to_vcf(g_ref_ancestrydna_b37, data_lines);
    //report_semaphore();
    return;
  }

  console.log("unknown format");
}

//                                         _              _ 
//  _ __ _ _ ___  __ ___ ______  _  _ _ __| |___  __ _ __| |
// | '_ \ '_/ _ \/ _/ -_|_-<_-< | || | '_ \ / _ \/ _` / _` |
// | .__/_| \___/\__\___/__/__/  \_,_| .__/_\___/\__,_\__,_|
// |_|                               |_|                    


function queue_example_genome() {
  start_spinner();
  setTimeout(function() {
    load_default_clinvar();
    load_example_vcf();
  }, 1);
}

function queue_clinvar(inp_id) {
  start_spinner();
  setTimeout(function() { process_clinvar_upload(inp_id); }, 1);
}


function process_clinvar_upload() {
  var inp_id = "clinvar_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(fn);

  g_clinvar_ready  = false;
  g_clinvar_queued = true;

  var reader = new FileReader();
  reader.onload = load_clinvar_file;
  reader.readAsDataURL(fn);

}

function process_ref_23andme_upload() {
  var inp_id = "ref_23andme_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(fn);

  g_ref_23andme_ready = false;

  var reader = new FileReader();
  reader.onload = load_ref_23andme_file;
  reader.readAsDataURL(fn);

}

function process_ref_ancestrydna_upload() {
  var inp_id = "ref_ancestrydna_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(fn);

  g_ref_ancestrydna_ready = false;

  var reader = new FileReader();
  reader.onload = load_ref_ancestrydna_file;
  reader.readAsDataURL(fn);

}

function queue_genome(inp_id) {

  console.log("... queueing genome");

  start_spinner();
  setTimeout(function() { process_genome_upload(inp_id); }, 1);
}

function process_genome_upload(inp_id) {
  //var inp_id = "genome_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(">>>>>>>>>>>>", fn);

  g_genome_ready = false;
  if (!g_clinvar_ready) {
    if (!g_clinvar_queued) {
      load_default_clinvar();
    }
  }

  var reader = new FileReader();
  reader.onload = load_genome_file;
  reader.readAsDataURL(fn);
}

function start_spinner() {
  $("#welcome_text").css("display", "none");
  $("#processing").css("display", "inline");
  $("#main").css("display", "none");
}

function remove_spinner() {
  $("#processing").css("display", "none");
  $("#main").css("display", "inline");
}

var g_first_sort = true;
var g_variants = [];

//  _      _ _   
// (_)_ _ (_) |_ 
// | | ' \| |  _|
// |_|_||_|_|\__|
//               


function init_table() {

  // https://mottie.github.io/tablesorter/docs/#Examples
  // lets try jquery tablesorter
  //
  $(function() {

    $.tablesorter.themes.bootstrap = {
			// these classes are added to the table. To see other table classes available,
			// look here: http://getbootstrap.com/css/#tables
			table        : 'table table-bordered table-striped',
			caption      : 'caption',
			// header class names
			header       : 'bootstrap-header', // give the header a gradient background (theme.bootstrap_2.css)
			sortNone     : '',
			sortAsc      : '',
			sortDesc     : '',
			active       : '', // applied when column is sorted
			hover        : '', // custom css required - a defined bootstrap style may not override other classes
			// icon class names
			icons        : '', // add "icon-white" to make them white; this icon class is added to the <i> in the header
			iconSortNone : 'bootstrap-icon-unsorted', // class name added to icon when column is not sorted
			iconSortAsc  : 'glyphicon glyphicon-chevron-up', // class name added to icon when column has ascending sort
			iconSortDesc : 'glyphicon glyphicon-chevron-down', // class name added to icon when column has descending sort
			filterRow    : '', // filter row class; use widgetOptions.filter_cssFilter for the input/select element
			footerRow    : '',
			footerCells  : '',
			even         : '', // even row zebra striping
			odd          : ''  // odd row zebra striping
    };


		// call the tablesorter plugin and apply the uitheme widget
		$("table").tablesorter({

			// this will apply the bootstrap theme if "uitheme" widget is included
			// the widgetOptions.uitheme is no longer required to be set
			theme : "bootstrap",
			widthFixed: true,
			headerTemplate : '{content} {icon}', // new in v2.7. Needed to add the bootstrap icon!

      sortList: [[2,0]],

      //headers: { 2 : { sortInitialOrder: 'asc' } },

			// widget code contained in the jquery.tablesorter.widgets.js file
			// use the zebra stripe widget if you plan on hiding any rows (filter widget)
			//widgets : [ "uitheme", "filter", "columns", "zebra" ],
			widgets : [ "uitheme", "filter", "zebra" ],
			widgetOptions : {
				// using the default zebra striping class name, so it actually isn't included in the theme variable above
				// this is ONLY needed for bootstrap theming if you are using the filter widget, because rows are hidden
				zebra : ["even", "odd"],

				// class names added to columns when sorted
				columns: [ "primary", "secondary", "tertiary" ],

				// reset filters button
				filter_reset : ".reset",

				// extra css class name (string or array) added to the filter element (input or select)
				filter_cssFilter: "form-control",

        filter_functions : {
          1 : {
            //"all" : function(e,n,f,i,$r,c,data) { return true; },
            "pathogenic" : function(e,n,f,i,$r,c,data) { return /^[pP]athogenic$/.test(e); },
            "probably pathogenic" : function(e,n,f,i,$r,c,data) { return /^[pP]robably pathogenic$/.test(e); },
            "affecting drug response" : function(e,n,f,i,$r,c,data) { return /^[aA]ffecting drug response$/.test(e); },
            "affecting histocompatibility" : function(e,n,f,i,$r,c,data) { return /^[aA]ffecting histocompatibility$/.test(e); },

            "non-pathogenic" : function(e,n,f,i,$r,c,data) { return /^[nN]on-pathogenic$/.test(e); },
            "probably non-pathogenic" : function(e,n,f,i,$r,c,data) { return /^[pP]robably non-pathogenic$/.test(e); },

            "unknown" : function(e,n,f,i,$r,c,data) { return /^[uU]nknown/.test(e); },
            "untested" : function(e,n,f,i,$r,c,data) { return /^[uU]tested/.test(e); }
          }
        }

				// set the uitheme widget to use the bootstrap theme class names
				// this is no longer required, if theme is set
				// ,uitheme : "bootstrap"

			}
		})
		.tablesorterPager({

			// target the pager markup - see the HTML block below
			container: $(".ts-pager"),

			// target the pager page select dropdown - choose a page
			cssGoto  : ".pagenum",

			// remove rows from the table to speed up the sort of large tables.
			// setting this to false, only hides the non-visible rows; needed if you plan to add/remove rows with the pager enabled.
			removeRows: false,

			// output string - default is '{page}/{totalPages}';
			// possible variables: {page}, {totalPages}, {filteredPages}, {startRow}, {endRow}, {filteredRows} and {totalRows}
			output: '{startRow} - {endRow} / {filteredRows} ({totalRows})'

		});

  });


}

//                     _ _           _   _          
//  __ ___  ___ _ _ __| (_)_ _  __ _| |_(_)___ _ _  
// / _/ _ \/ _ \ '_/ _` | | ' \/ _` |  _| / _ \ ' \ 
// \__\___/\___/_| \__,_|_|_||_\__,_|\__|_\___/_||_|
//                                                  


function report_semaphore() {

  // Make sure our genome is ready and clinvar data is
  // ready.
  //
  if ((!g_genome_ready) || (!g_clinvar_ready)) {
    return;
  }

  // Clinvar significance code mapping
  //
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

  // Run the report
  //
  var data = vcf2clinvar.clinvar_report_json(g_clinvar_lines, g_genome_lines);
  var res = data.results;

  // Remove rows in the table.
  // Be careful, more work needs to be done,
  // see below
  //
  var chld = $("#table_body").children();
  for (var ii=0; ii<chld.length; ii++) {
    chld[ii].remove();
  }
  
  // I think jquery or whatever else we use (tablesorter?) has a cache
  // that doesn't update until you tell it explicitely to.
  //
  // Clear the cache.
  //
  $("#table_body").trigger("addRows").trigger("update").trigger("appendCache").trigger("applyWidgets");

  var rows = [];
  for (var ii=0; ii<res.length; ii++) {
    var s =  [ res[ii][4], _sig_lookup[res[ii][5]], res[ii][6], res[ii][7] ].join("</td><td>") ;
    rows.push(s);
  }

  var row_str = "<tr><td>" + rows.join("</td></tr><tr><td>") + "</td></tr>";

  var $row = $(row_str);

  var body = $("#table_body");
  body.append($row).trigger("addRows", [$row, false]);

  // Do the default filter and sorts
  //
  $("#table").find("select").val("pathogenic");
  $("#table").find("th:contains(Significance)").trigger("search");

  // This needs work...
  // We need to do this more intelligently:
  // https://mottie.github.io/tablesorter/docs/example-trigger-sort.html
  // https://mottie.github.io/tablesorter/docs/example-triggers.html
  //
  //$("#table").find("th:contains(Allele)").trigger("sort");
  //
  // oof: http://wowmotty.blogspot.com/2011/06/jquery-tablesorter-missing-docs.html
  //
  $("#table").trigger("sorton", [[[2,0]]] );

  remove_spinner();
}

function data_23andme_semaphore() {

  console.log("data_23andme_semaphore");

  if (g_process_data) {

    console.log("data_23andme_semaphore ... still needs processing");

    if ((typeof g_data_lines !== "undefined") && (g_ref_23andme_ready)) {

      console.log("data_23andme_semaphore ... processing");

      g_genome_lines = convert_23andme_to_vcf(g_ref_23andme_b37, g_data_lines);
      g_genome_ready = true;

      report_semaphore();

      g_process_data = false;
    }

    else if (!g_ref_23andme_ready) { load_default_ref_23andme(); }

  }

}

function data_ancestrydna_semaphore() {

  console.log("data_ancestrydna_semaphore");

  if (g_process_data) {

    console.log("data_ancestrydna_semaphore ... still needs processing");

    if ((typeof g_data_lines !== "undefined") && (g_ref_ancestrydna_ready)) {

      console.log("data_ancestrydna_semaphore ... processing");

      g_genome_lines = convert_ancestrydna_to_vcf(g_ref_ancestrydna_b37, g_data_lines);
      g_genome_ready = true;

      report_semaphore();

      g_process_data = false;

    }

    else if (!g_ref_ancestrydna_ready) { load_default_ref_ancestrydna(); }

  }

}



//     _      __           _ _     _              _ 
//  __| |___ / _|__ _ _  _| | |_  | |___  __ _ __| |
// / _` / -_)  _/ _` | || | |  _| | / _ \/ _` / _` |
// \__,_\___|_| \__,_|\_,_|_|\__| |_\___/\__,_\__,_|
//                                                  


function load_default_clinvar() {
  g_clinvar_ready = false;
  g_clinvar_queued = true;

  console.log("loading default clinvar");

  // Load the ClinVar (stripped down and compressed) VCF
  // as a background process.
  //
  var cln_wurk = new Worker("js/clinvar-worker.js");
  cln_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    g_clinvar_lines = raw_data.split("\n");

    g_clinvar_ready = true;
    g_clinvar_queued = false;

    report_semaphore();
    console.log("clinvar:", g_clinvar_lines.length, "lines loaded");
  });
}

function load_default_ref_23andme() {

  g_ref_23andme_ready = false;
  g_ref_23andme_queued = true;

  // Load the reference positions for 23andMe as a
  // background process.
  //
  var ttam_ref_wurk = new Worker("js/ref-b37-worker.js");
  ttam_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_ref_23andme_b37, raw_lines);

    g_ref_23andme_ready = true;
    g_ref_23andme_queued = false;

    console.log("g_ref_23andme loaded");
    data_23andme_semaphore();
  });
  ttam_ref_wurk.postMessage("../data/23andme_reference_b37.txt.gz");
}

function load_default_ref_ancestrydna() {

  g_ref_ancestrydna_ready = false;
  g_ref_ancestrydna_queued = true;


  // Load the reference positions for AncestryDNA as a
  // background process.
  //
  var adna_ref_wurk = new Worker("js/ref-b37-worker.js");
  adna_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_ref_ancestrydna_b37, raw_lines);

    g_ref_ancestrydna_ready = true;
    g_ref_ancestrydna_queued = false;

    console.log("g_ref_ancestrydna loaded");

    data_ancestrydna_semaphore();
  });
  adna_ref_wurk.postMessage("../data/ancestrydna_reference_b37.txt.gz");
}

function load_example_vcf() {

  g_genome_ready = false;

  // Load the example VCF as a background process
  //
  var example_wurk = new Worker("js/example-vcf-worker.js");
  example_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    g_genome_lines = raw_data.split("\n");

    g_genome_ready = true;

    report_semaphore();

    console.log("example:", g_genome_lines.length, "lines loaded");
  });
}


