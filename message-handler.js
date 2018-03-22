//------------------------------------------------------------------------------------------------------------------------
function setupMessageHandlers()
{
    Shiny.addCustomMessageHandler("roi",
                                  function(message) {
                                      console.log("roi message");
                                      console.log(message);
                                      window.igvBrowser.search(message.roi);
                                  }
                                 );

    Shiny.addCustomMessageHandler("displaySnps",

                                  function(message){
                                      console.log("displaySnps message");
                                      console.log(message);
                                      var trackName = "SNPs";
                                      var bedFileName = message.filename;
                                      var displayMode = "EXPANDED";
                                      var color = "red";
                                      var minValue = 0;
                                      var maxValue = 10;
                                      var trackHeight = 120;
                                      var url = window.location.href + bedFileName;
                                      var config = {format: "bedgraph",
                                                    name: trackName,
                                                    url: url,
                                                    min: minValue,
                                                    max: maxValue,
                                                    indexed: false,
                                                    displayMode: displayMode,
                                                    sourceType: "file",
                                                    color: color,
                                                    height: trackHeight,
                                                    type: "wig"};
                                      console.log("=== about to loadTrack");
                                      console.log(config)
                                      window.igvBrowser.loadTrack(config);
                                  })

    Shiny.addCustomMessageHandler("displayBedTrack", displayBedTrack);

} // setupMessageHandlers
//------------------------------------------------------------------------------------------------------------------------
function displayBedTrack(message)
{
  console.log("displayBedTrack, message:");
  console.log(message);

  var url = window.location.href + message.filename;

  var config = {format: "bed",
                name: message.trackName,
                url: url,
                indexed: false,
                displayMode: "SQUISHED",
                sourceType: "file",
                color: message.color,
                height: message.trackHeight,
                type: "annotation"};
  console.log("=== about to loadTrack");
  console.log(config)
  window.igvBrowser.loadTrack(config);



} // displayBedTrack
//------------------------------------------------------------------------------------------------------------------------
setTimeout(setupMessageHandlers, 5000);

