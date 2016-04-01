//Copyright (c) 2016 eXact Lab s.r.l.
//Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

//Export functions. Change these to modify how the polygons are saved.

//Selects the content of export div
function selectText(e) {

    var node = document.getElementById('export');

    if ( document.selection ) {
        var range = document.body.createTextRange();
        range.moveToElementText( node  );
        range.select();
    } else if ( window.getSelection ) {
        var range = document.createRange();
        range.selectNodeContents( node );
        window.getSelection().removeAllRanges();
        window.getSelection().addRange( range );
    }
}

//Exports a point list into a snippet of code for Bit.Sea
function exportPointsList(points) {
    $('#export').empty();
    $('#ui').empty();
    $('#ui').append("<input type='button' value='Select text' id='expSelAllbtn'></input>\n");
    $('#ui').append("<input type='button' value='Save all' id='saveAllbtn'></input>\n");
    $('#expSelAllbtn').click(selectText);
    $('#saveAllbtn').click(saveAllPoly);
    $('#export').append('p = Polygon([');
    var lat_list = [];
    var lon_list = [];
    for (i in points) {
        lat_list.push(points[i].lat);
        lon_list.push(points[i].lng);
    }
    $('#export').append(lon_list.toString());
    $('#export').append('],[');
    $('#export').append(lat_list.toString());
    $('#export').append('])');
}

//Save a text file with all the polygons
function saveAllPoly(e) {
    var fullText = "";
    polys = drawnItems.getLayers();
    for (i in polys) {
        p = polys[i];
        points = p.getLatLngs();
        fullText = fullText + "p" + i + " = Polygon([";
        var lat_list = [];
        var lon_list = [];
        for (i in points) {
            lat_list.push(points[i].lat);
            lon_list.push(points[i].lng);
        }

        fullText = fullText + lon_list.toString() + "],["
        fullText = fullText + lat_list.toString() + "])\n"
    }
    var blob = new Blob([fullText], {type: "text/plain;charset=utf-8"});
    saveAs(blob, "polygon.py");
}
