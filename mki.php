<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1,maximum-scale=1" />
        <link rel="stylesheet" href="styles.css">
        <title>MarkInChi experiments</title>
    </head>
    <body>


        <section class="navbar">
            <h1 class="title hidden">Navigation bar</h1>
            <a href="./index.html" class="active">MarkInChI</a>
            <div class="dropdown">
                <button class="dropdownButton">InChI Resources &#9660;</button>
                <div class="dropdownContent">
                    <a href="https://www.inchi-trust.org">InChI Trust</a>
		    <a href="https://iupac.org/who-we-are/committees/committee-details/?body_code=802">IUPAC InChI subcommittee</a>
                    <a href="http://www-rinchi.ch.cam.ac.uk/qrinchi/">QRInChI</a>
                    <a href="https://doi.org/10.1186/s13321-018-0277-8">RInChI: the reaction InChI</a>
                </div>
            </div>
            <div class="dropdown">
                <button class="dropdownButton">About this program &#9660;</button>
                <div class="dropdownContent">
                    <a href="http://orcid.org/0000-0002-8693-9136">J. M. Goodman</a>
                    <a href="http://orcid.org/0000-0002-2824-1805">E. D. J. Goodman</a>
                    <a href="https://github.com/EdmundGoodman/InChiMarkush">GitHub Source Code</a>
                </div>
            </div>
        </section>



        <section id="contentWrapper">
            <h1 class="title">
                Program output
            </h1>

            <h2 class="heading">
		The deconstructed MarkInChi output is:
        </h2>

        <textarea id="output", rows="5">

<?php
if ($_GET['InChIKey']) {
    $code = trim($_GET['InChIKey'], " ");

    //Regex schema for the InChi key - reports errors helpfully, and also prevents
    //security issues by severly restricting the valid input domain (need ^ and $
    //to anchor to the entire string)
    $InChiSchema = "/^InChI=1[a-zA-Z0-9\/\>\<\!\,\-\(\)]*$/";

    if (preg_match($InChiSchema, $code)) {
        $escapedCode = escapeshellarg($code);

        //If this isn't the raw posted code plus quotation marks,
        //something malicious is going on, so redirect to the error page
        $expectedEscapedCode = "'" . $code . "'";
        if ($expectedEscapedCode != $escapedCode) {
            header("Location: ./error.html");
            exit();
        }

        //Run the script to deconstruct the MarkInChi
        $safeCommand = "python3 markinchi2inchi.py " . $escapedCode;
        passthru($safeCommand);

    } else {
        //The code is of an invalid schema, so redirect to the error page
	header("Location: ./error.html");
        exit();
    }

} else {
	//The code field isn't set, so redirect to the error page
    	header("Location: ./error.html");
    	exit();
}
?>


        </textarea>
        <br>
        <button id="copybutton" onclick="myFunction()" onmouseout="outFunc()">
            Copy output to clipboard
        </button>
        <script>
        function myFunction() {
            document.querySelector("textarea").select();
            document.execCommand('copy');
        }
        </script>


        
        </section>

    </body>
</html>
