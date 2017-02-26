importScripts('gunzip.min.js');

var xhr = new XMLHttpRequest();
xhr.open('GET', '../data/clinvar.vcf.gz', true);
xhr.responseType = 'arraybuffer';

var g_self = self;

xhr.onload = function(e) {

  console.log("clinvar ok>>");

  var uInt8Array = new Uint8Array(this.response);
  var gunzip = new Zlib.Gunzip(uInt8Array);
  var unpacked_uInt8Array = gunzip.decompress();
  g_self.postMessage(unpacked_uInt8Array);
};

xhr.send();


