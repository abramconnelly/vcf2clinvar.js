#!/usr/bin/env node

var g_config = {
  "genomeBuild" : "b37",
  "outputType" : "csv"
};

var g_genomeBuild = "b37";
var g_fn = undefined;

var fs = require("fs"),
    program = require("commander"),
    vcf2clinvar = require("./js/vcf2clinvar.js");

var indent = "\n\t\t\t\t\t"

program
  .arguments('<file>')
  .option('-c, --clinvarfile <clinvarfile>', 'ClinVar VCF file')
  /*
  .option('-C, --clinvardir <clinvardir>', 'ClinVar VCF directory (either this or -c must be' +
                                           'specified). This option will use' +
                                           'vcf2clinvar.clinvar_update to automatically check and' +
                                           'import the most recent ClinVar file to this directory.')
                                           */
  //.option('-i, --input <input>', "Input VCF file ['.vcf', '.vcf.gz', '.vcf.bz2']." +
  //                               'Uncompressed genome data is also accepted via stdin.')
  .option('-i, --input <input>', "Input VCF file")
  .option('-t, --type <type>', "Output report type ('csv' or 'json'). Defaults to csv." +
                               indent + "CSV Report: Reports all genome variants matching" +
                               indent + "ClinVar records, and some summary ClinVar data from" +
                               indent + "these records. Header lines with metadata begin with" +
                               indent + "'##'. JSON Report: Reports genome variants matching" +
                               indent + "ClinVar records (no record information is included).")
  .option('-n, --notes <notes>', 'Notes (JSON format) to include in report. (JSON report only)')
  .option('-g, --genome-build <genomebuild>', "Genome build to include in report ('b37' or 'b38').")
  .action(function(file) {
    console.log('clinvarfile %s, input %s, type %s, notes %s, genome-build %s, file %s',
        program.clinvarfile, program.input, program.type, program.notes, program.genomeBuild, file);
    g_config["input"] = file;
  })
  .parse(process.argv)

if (typeof program.clinvarfile === "undefined") {
  console.log("provide clinvarfile");
  program.help();
}

g_config["clinvarfile"] = program.clinvarfile;
g_config["notes"] = program.notes;

if (typeof program.type !== "undefined") {
  g_config["outputType"] = program.type;
}

if ((g_config["outputType"] != "csv") && (g_config["outputType"] != "json")) {
  console.log("output type must be 'csv' or 'json'");
  program.help();
}

if ((typeof program.input === "undefined") && (typeof g_config["input"] === "undefined")) {
  console.log("specify only one input VCF file");
  program.help();
}

if (typeof program.input !== "undeifned") {
  g_config["input"] = program.input;
}

if (typeof g_config["input"] === "undefined" ) {
  console.log("specify input VCF file");
  program.help();
}

if (typeof program.genomeBuild !== "undefined") {
  g_config["genomeBuild"] = program.genomeBuild;
}

if ((g_config["genomeBuild"] != "b37") && (g_config["genomeBuild"] != "b38")) {
  console.log("genome build must be b37 or b38 (defaults to b37)");
  program.help();
}


console.log(g_config);

console.log('... clinvarfile %s, input %s, type %s, notes %s, genome-build %s, file %s',
    program.clinvarfile, program.input, program.type, program.notes, program.genomeBuild, g_fn)

var clinvar_lines = fs.readFileSync(g_config.clinvarfile, {encoding:"utf8"}).split("\n");
var input_lines = fs.readFileSync(g_config.input, {encoding:"utf8"}).split("\n");

var res = vcf2clinvar.clinvar_report(clinvar_lines, input_lines);
console.log(res.join("\n"));
