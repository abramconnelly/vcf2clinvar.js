importScripts('gunzip.min.js');

var g_self = self;
var xhr = new XMLHttpRequest();

xhr.onload = function(e) {
  var uInt8Array = new Uint8Array(this.response);
  var gunzip = new Zlib.Gunzip(uInt8Array);
  var unpacked_uInt8Array = gunzip.decompress();
  g_self.postMessage(unpacked_uInt8Array);
};

self.addEventListener("message", function(e) {

  console.log("ref-b37-worker: got", e.data);

  //xhr.open('GET', '../data/23andme_reference_b37.txt.gz', true);
  xhr.open('GET', e.data, true);
  xhr.responseType = 'arraybuffer';
  xhr.send();
});



