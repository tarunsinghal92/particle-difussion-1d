<?php

/**
 * url
 */
define('SITEURL', 'particle-difussion-1d/');

/**
 * time step (seconds)
 * sample value : 1 day = 1 * 8.64 * 10^4 seconds
 */
define('TIME_STEP', (1 * 86400));

/**
 * initial Concenteration (in percent)
 */
define('INITIAL_CONC', 1.0);

/**
 * overall time (seconds)
 * sample value : 100 day
 */
define('OVERALL_TIME', (100 * 86400));

/**
 * deltax of material (m)
 */
define('STEP_SIZE', 0.02);

/**
 * total depth 100m or 0.1m
 */
define('TOTAL_LENGTH', 0.1);

/**
 * diffusion constant (m2/s)
 */
define('DIFFUSION_CONSTANT', (0.000000000006));

/**
 * theta (Crank nicolson factor)
 */
define('THETA', (0.5));

?>
