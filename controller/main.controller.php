<?php

/**
 * Main Controller Class
 */

class MainController extends Common{

    //all results
    private $results;

    //default function
    public function __construct(){

        //loop for various values
        $dt = [0.1, 0.2, 0.5, 1.0, 5.0, 10, 20];
        $dx = [0.001, 0.002, 0.005, 0.01, 0.02]; 

        echo '#,' . implode(',', $dt) . '<br>';
        foreach ($dx as $k2 => $xval) {
            echo $xval ;
            foreach ($dt as $k1 => $tval) {
                //get the results
                $test = new DiffusionModel($xval, ($tval * 86400)); // detaX in m, deltaT in sec
                $test->run();

                //results
                echo ',' . $test->error_results;
                ob_flush();flush();
            }
            echo '<br>';
        }
        die;
    }

    public function show_template(){

        //show template
        require_once 'views/main.phtml';
    }
}

?>
