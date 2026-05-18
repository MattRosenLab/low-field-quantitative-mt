function [M]=Read_Tecmag_hdr()
M=['TNT1000                    ',' 000000 ',' char   ',' 008 ','Version ID                           ';
   'TMAG                       ',' 000008 ',' char   ',' 004 ','Tecmag Tag                           ';
   'BOOL                       ',' 000012 ',' short  ',' 001 ','BOOlean Value                        ';
   'LEN                        ',' 000016 ',' short  ',' 001 ','Length of Tecmag structure           ';

% Number of points and scans in all dimensions
   'npts                       ',' 000020 ',' long   ',' 004 ','                                     ';
   'actual_npts                ',' 000036 ',' long   ',' 004 ','                                     '; 
   'acq_points                 ',' 000052 ',' long   ',' 001 ','                                     ';
   'npts_start                 ',' 000056 ',' long   ',' 004 ','                                     ';
   'scans                      ',' 000068 ',' long   ',' 001 ','                                     ';
   'actual_scans               ',' 000072 ',' long   ',' 001 ','                                     ';
   'dummy_scans                ',' 000076 ',' long   ',' 001 ','                                     ';
   'repeat_times               ',' 000080 ',' long   ',' 001 ','                                     ';
   'sadimension                ',' 000084 ',' long   ',' 001 ','                                     ';
%   'space1                     ',' 000088 ',' char   ',' 004 ','                                     ';

% Field and frequencies
   'magnet_field               ',' 000096 ',' float64',' 001 ','                                     ';
   'ob_freq                    ',' 000104 ',' float64',' 004 ','                                     ';
   'base_freq                  ',' 000136 ',' float64',' 004 ','                                     ';
   'offset_freq                ',' 000168 ',' float64',' 004 ','                                     ';
   'ref_freq                   ',' 000200 ',' float64',' 001 ','                                     ';
   'NMR_freq                   ',' 000208 ',' float64',' 001 ','                                     ';
%   'spaces2                    ',' 000216 ',' char   ',' 044 ','                                     ';
%
% spectral width, dwell and filter
   'sw                         ',' 000260 ',' float64',' 004 ','                                     ';
   'dwell                      ',' 000292 ',' float64',' 004 ','                                     ';
   'filter                     ',' 000324 ',' float64',' 001 ','                                     ';
   'experiment_time            ',' 000332 ',' float64',' 001 ','                                     ';
   'acq_time                   ',' 000340 ',' float64',' 001 ','                                     ';
%
   'last_delay                 ',' 000348 ',' float64',' 001 ','                                     ';
%
   'spectrum_direction         ',' 000356 ',' short  ',' 001 ','                                     ';
   'hardware_sideband          ',' 000358 ',' short  ',' 001 ','                                     ';
   'Taps                       ',' 000360 ',' short  ',' 001 ','                                     ';
   'Type                       ',' 000362 ',' short  ',' 001 ','                                     ';
   'dDigRec                    ',' 000364 ',' ubit4  ',' 001 ','                                     ';
   'nDigitalCenter             ',' 000368 ',' long   ',' 001 ','                                     ';
%   'spaces3                    ',' 000372 ',' char   ',' 016 ','                                     ';
%
% Hardware settings
   'tramsmitter_gain           ',' 000388 ',' short  ',' 001 ','                                     ';
   'receiver_gain              ',' 000390 ',' short  ',' 001 ','                                     ';
%   'spaces4                    ',' 000392 ',' char   ',' 016 ','                                     ';
%
% spinning speed information
   'set_spin_rate              ',' 000408 ',' ushort ',' 001 ','                                     ';
   'actual_spin_rate           ',' 000410 ',' ushort ',' 001 ','                                     ';
%
% Lock information
   'lock_field                 ',' 000412 ',' short  ',' 001 ','                                     ';
   'lock_power                 ',' 000414 ',' short  ',' 001 ','                                     ';
   'lock_gain                  ',' 000416 ',' short  ',' 001 ','                                     ';
   'lock_phase                 ',' 000418 ',' short  ',' 001 ','                                     ';
   'lock_freq_mhz              ',' 000420 ',' float64',' 001 ','                                     ';
   'lock_ppm                   ',' 000428 ',' float64',' 001 ','                                     ';
   'H2O_freq_ref               ',' 000436 ',' float64',' 001 ','                                     ';
%   'spaces5                    ',' 000344 ',' char   ',' 016 ','                                     ';
%
% VT information
   'set_temperature            ',' 000460 ',' float64',' 001 ','                                     ';
   'actual_temperature         ',' 000468 ',' float64',' 001 ','                                     ';
%
% Shim information
   'shim_unit                  ',' 000476 ',' float64',' 001 ','                                     ';
   'shims                      ',' 000484 ',' short  ',' 036 ','                                     ';
   'shims_FWHM                 ',' 000556 ',' float64',' 001 ','                                     ';
%
% Hardware specific information
   'HH_dcpl_attn               ',' 000564 ',' short  ',' 001 ','                                     ';
%
   'DF_DN                      ',' 000566 ',' short  ',' 001 ','                                     ';
   'F1_tran_mode               ',' 000568 ',' short  ',' 007 ','                                     ';
   'dec_BW                     ',' 000582 ',' short  ',' 001 ','                                     ';
%
   'grd_orientation            ',' 000584 ',' char   ',' 004 ','                                     ';
%   'spaces6                    ',' 000588 ',' char   ',' 296 ','                                     ';
%
% text variables
   'date                       ',' 000884 ',' char   ',' 032 ','                                     ';
   'nucleus                    ',' 000916 ',' char   ',' 016 ','                                     ';
   'nucleus_2D                 ',' 000932 ',' char   ',' 016 ','                                     ';
   'nucleus_3D                 ',' 000948 ',' char   ',' 016 ','                                     ';
   'nucleus_4D                 ',' 000964 ',' char   ',' 016 ','                                     ';
   'sequence                   ',' 000980 ',' char   ',' 032 ','                                     ';
   'lock_solvent               ',' 001012 ',' char   ',' 016 ','                                     ';
   'lock_nucleus               ',' 001028 ',' char   ',' 016 ','                                     '];
