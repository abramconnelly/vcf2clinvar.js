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

function load_file(e) {
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

  g_genome_lines = undefined;

  var data_lines = data_str.split("\n");
  if (data_lines.length==0) { return; }
  if (data_lines[0].search(/##fileformat=VCF/)>=0) {
    g_genome_lines = data_lines;
    report_semaphore();
    return;
  }

  if (data_lines[0].search(/23and[mM]e/)>=0) {
    console.log(">>> 23andme data?");

    //var vcf_lines = convert_23andme_to_vcf(g_23andme_ref_b37, data_lines);
    g_genome_lines = convert_23andme_to_vcf(g_23andme_ref_b37, data_lines);
    report_semaphore();
    return;
  }

  if (data_lines[0].search(/Ancestry[dD][nN][aA]/)>=0) {
    console.log(">>> Ancestry DNA?");

    //var vcf_lines = convert_ancestrydna_to_vcf(g_ancestrydna_ref_b37, data_lines);
    g_genome_lines = convert_ancestrydna_to_vcf(g_ancestrydna_ref_b37, data_lines);
    report_semaphore();
    return;
  }

  console.log("unknown format");
}

function process_genome_upload() {
  var inp_id = "genome_upload";
  var fn = document.getElementById(inp_id).files[0];
  console.log(fn);

  var reader = new FileReader();
  reader.onload = load_file;
  reader.readAsDataURL(fn);
}

function remove_spinner() {
  $("#spinner").css("display", "none");
  $("#main").css("display", "inline");

}

var g_first_sort = true;
var g_variants = [];

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

function on_ready() {
  init_table();

  //TODO: load this staggered or only on user initiation
  //

  // Load the ClinVar (stripped down and compressed) VCF
  // as a background process.
  //
  var cln_wurk = new Worker("js/clinvar-worker.js");
  cln_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    g_clinvar_lines = raw_data.split("\n");

    g_clinvar_ready = true;

    report_semaphore();
    console.log("clinvar:", g_clinvar_lines.length, "lines loaded");
  });

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

  // Load the reference positions for 23andMe as a
  // background process.
  //
  var ttam_ref_wurk = new Worker("js/ref-b37-worker.js");
  ttam_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_23andme_ref_b37, raw_lines);
    console.log("g_23andme_ref loaded");
  });
  ttam_ref_wurk.postMessage("../data/23andme_reference_b37.txt.gz");

  // Load the reference positions for AncestryDNA as a
  // background process.
  //
  var adna_ref_wurk = new Worker("js/ref-b37-worker.js");
  adna_ref_wurk.addEventListener("message", function(e) {
    var raw_data = new TextDecoder("utf-8").decode(e.data);
    var raw_lines = raw_data.split("\n");
    load_gt_reference(g_ancestrydna_ref_b37, raw_lines);
    console.log("g_ancestrydna_ref loaded");
  });
  adna_ref_wurk.postMessage("../data/ancestrydna_reference_b37.txt.gz");

}
