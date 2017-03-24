/*
  File autogenerated by gengetopt version 2.22.6
  generated with the following command:
  gengetopt --file-name=cmdline 

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.h"

const char *gengetopt_args_info_purpose = "";

const char *gengetopt_args_info_usage = "Usage: zeno [OPTIONS]...";

const char *gengetopt_args_info_versiontext = "";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_full_help[] = {
  "  -h, --help                    Print help and exit",
  "      --full-help               Print help, including hidden options, and exit",
  "  -V, --version                 Print version and exit",
  "  -i, --input-file=STRING       Input file name",
  "      --num-walks=LONGLONG      Number of walk-on-spheres walks to perform",
  "      --num-interior-samples=LONGLONG\n                                Number of interior samples to take",
  "      --max-rsd-capacitance=DOUBLE\n                                Perform walk-on-spheres walks until the\n                                  relative standard deviation of the\n                                  capacitance drops below this value.  Relative\n                                  standard deviation is defined as\n                                  (Standard_Deviation/Mean)*100%",
  "      --max-rsd-polarizability=DOUBLE\n                                Perform walk-on-spheres walks until the\n                                  relative standard deviation of the mean\n                                  electric polarizability drops below this\n                                  value.  Relative standard deviation is\n                                  defined as (Standard_Deviation/Mean)*100%",
  "      --max-rsd-volume=DOUBLE   Take interior samples until the relative\n                                  standard deviation of volume drops below this\n                                  value.  Relative standard deviation is\n                                  defined as (Standard_Deviation/Mean)*100%",
  "      --min-num-walks=LONGLONG  Minimum number of walk-on-spheres walks to\n                                  perform when using max-rsd stopping\n                                  conditions  (default=`1000')",
  "      --min-num-interior-samples=LONGLONG\n                                Minimum number of interior samples to take when\n                                  using max-rsd stopping conditions\n                                  (default=`10000')",
  "      --compute-form            Compute form factor",
  "      --num-threads=INT         Number of threads to use  (default=Number of\n                                  logical cores)",
  "      --seed=INT                Seed for the random number generator\n                                  (default=Randomly set)",
  "      --frac-error-bound=DOUBLE Fractional error bound for nearest neighbor\n                                  search  (default=`0')",
  "      --surface-points-file=STRING\n                                Name of file for writing the surface points\n                                  from Walk-on-Spheres",
  "      --interior-points-file=STRING\n                                Name of file for writing the interior sample\n                                  points",
  "      --print-counts            Print statistics related to counts of hit\n                                  points",
  "      --print-benchmarks        Print detailed RAM and timing information",
    0
};

static void
init_help_array(void)
{
  gengetopt_args_info_help[0] = gengetopt_args_info_full_help[0];
  gengetopt_args_info_help[1] = gengetopt_args_info_full_help[1];
  gengetopt_args_info_help[2] = gengetopt_args_info_full_help[2];
  gengetopt_args_info_help[3] = gengetopt_args_info_full_help[3];
  gengetopt_args_info_help[4] = gengetopt_args_info_full_help[4];
  gengetopt_args_info_help[5] = gengetopt_args_info_full_help[5];
  gengetopt_args_info_help[6] = gengetopt_args_info_full_help[6];
  gengetopt_args_info_help[7] = gengetopt_args_info_full_help[7];
  gengetopt_args_info_help[8] = gengetopt_args_info_full_help[8];
  gengetopt_args_info_help[9] = gengetopt_args_info_full_help[9];
  gengetopt_args_info_help[10] = gengetopt_args_info_full_help[10];
  gengetopt_args_info_help[11] = gengetopt_args_info_full_help[11];
  gengetopt_args_info_help[12] = gengetopt_args_info_full_help[12];
  gengetopt_args_info_help[13] = gengetopt_args_info_full_help[13];
  gengetopt_args_info_help[14] = gengetopt_args_info_full_help[15];
  gengetopt_args_info_help[15] = gengetopt_args_info_full_help[16];
  gengetopt_args_info_help[16] = gengetopt_args_info_full_help[17];
  gengetopt_args_info_help[17] = gengetopt_args_info_full_help[18];
  gengetopt_args_info_help[18] = 0; 
  
}

const char *gengetopt_args_info_help[19];

typedef enum {ARG_NO
  , ARG_STRING
  , ARG_INT
  , ARG_DOUBLE
  , ARG_LONGLONG
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);

static int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error);

static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->full_help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->input_file_given = 0 ;
  args_info->num_walks_given = 0 ;
  args_info->num_interior_samples_given = 0 ;
  args_info->max_rsd_capacitance_given = 0 ;
  args_info->max_rsd_polarizability_given = 0 ;
  args_info->max_rsd_volume_given = 0 ;
  args_info->min_num_walks_given = 0 ;
  args_info->min_num_interior_samples_given = 0 ;
  args_info->compute_form_given = 0 ;
  args_info->num_threads_given = 0 ;
  args_info->seed_given = 0 ;
  args_info->frac_error_bound_given = 0 ;
  args_info->surface_points_file_given = 0 ;
  args_info->interior_points_file_given = 0 ;
  args_info->print_counts_given = 0 ;
  args_info->print_benchmarks_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->input_file_arg = NULL;
  args_info->input_file_orig = NULL;
  args_info->num_walks_orig = NULL;
  args_info->num_interior_samples_orig = NULL;
  args_info->max_rsd_capacitance_orig = NULL;
  args_info->max_rsd_polarizability_orig = NULL;
  args_info->max_rsd_volume_orig = NULL;
  args_info->min_num_walks_arg = 1000;
  args_info->min_num_walks_orig = NULL;
  args_info->min_num_interior_samples_arg = 10000;
  args_info->min_num_interior_samples_orig = NULL;
  args_info->num_threads_orig = NULL;
  args_info->seed_orig = NULL;
  args_info->frac_error_bound_arg = 0;
  args_info->frac_error_bound_orig = NULL;
  args_info->surface_points_file_arg = NULL;
  args_info->surface_points_file_orig = NULL;
  args_info->interior_points_file_arg = NULL;
  args_info->interior_points_file_orig = NULL;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{

  init_help_array(); 
  args_info->help_help = gengetopt_args_info_full_help[0] ;
  args_info->full_help_help = gengetopt_args_info_full_help[1] ;
  args_info->version_help = gengetopt_args_info_full_help[2] ;
  args_info->input_file_help = gengetopt_args_info_full_help[3] ;
  args_info->num_walks_help = gengetopt_args_info_full_help[4] ;
  args_info->num_interior_samples_help = gengetopt_args_info_full_help[5] ;
  args_info->max_rsd_capacitance_help = gengetopt_args_info_full_help[6] ;
  args_info->max_rsd_polarizability_help = gengetopt_args_info_full_help[7] ;
  args_info->max_rsd_volume_help = gengetopt_args_info_full_help[8] ;
  args_info->min_num_walks_help = gengetopt_args_info_full_help[9] ;
  args_info->min_num_interior_samples_help = gengetopt_args_info_full_help[10] ;
  args_info->compute_form_help = gengetopt_args_info_full_help[11] ;
  args_info->num_threads_help = gengetopt_args_info_full_help[12] ;
  args_info->seed_help = gengetopt_args_info_full_help[13] ;
  args_info->frac_error_bound_help = gengetopt_args_info_full_help[14] ;
  args_info->surface_points_file_help = gengetopt_args_info_full_help[15] ;
  args_info->interior_points_file_help = gengetopt_args_info_full_help[16] ;
  args_info->print_counts_help = gengetopt_args_info_full_help[17] ;
  args_info->print_benchmarks_help = gengetopt_args_info_full_help[18] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);

  if (strlen(gengetopt_args_info_versiontext) > 0)
    printf("\n%s\n", gengetopt_args_info_versiontext);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_print_full_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_full_help[i])
    printf("%s\n", gengetopt_args_info_full_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{

  free_string_field (&(args_info->input_file_arg));
  free_string_field (&(args_info->input_file_orig));
  free_string_field (&(args_info->num_walks_orig));
  free_string_field (&(args_info->num_interior_samples_orig));
  free_string_field (&(args_info->max_rsd_capacitance_orig));
  free_string_field (&(args_info->max_rsd_polarizability_orig));
  free_string_field (&(args_info->max_rsd_volume_orig));
  free_string_field (&(args_info->min_num_walks_orig));
  free_string_field (&(args_info->min_num_interior_samples_orig));
  free_string_field (&(args_info->num_threads_orig));
  free_string_field (&(args_info->seed_orig));
  free_string_field (&(args_info->frac_error_bound_orig));
  free_string_field (&(args_info->surface_points_file_arg));
  free_string_field (&(args_info->surface_points_file_orig));
  free_string_field (&(args_info->interior_points_file_arg));
  free_string_field (&(args_info->interior_points_file_orig));
  
  

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->full_help_given)
    write_into_file(outfile, "full-help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->input_file_given)
    write_into_file(outfile, "input-file", args_info->input_file_orig, 0);
  if (args_info->num_walks_given)
    write_into_file(outfile, "num-walks", args_info->num_walks_orig, 0);
  if (args_info->num_interior_samples_given)
    write_into_file(outfile, "num-interior-samples", args_info->num_interior_samples_orig, 0);
  if (args_info->max_rsd_capacitance_given)
    write_into_file(outfile, "max-rsd-capacitance", args_info->max_rsd_capacitance_orig, 0);
  if (args_info->max_rsd_polarizability_given)
    write_into_file(outfile, "max-rsd-polarizability", args_info->max_rsd_polarizability_orig, 0);
  if (args_info->max_rsd_volume_given)
    write_into_file(outfile, "max-rsd-volume", args_info->max_rsd_volume_orig, 0);
  if (args_info->min_num_walks_given)
    write_into_file(outfile, "min-num-walks", args_info->min_num_walks_orig, 0);
  if (args_info->min_num_interior_samples_given)
    write_into_file(outfile, "min-num-interior-samples", args_info->min_num_interior_samples_orig, 0);
  if (args_info->compute_form_given)
    write_into_file(outfile, "compute-form", 0, 0 );
  if (args_info->num_threads_given)
    write_into_file(outfile, "num-threads", args_info->num_threads_orig, 0);
  if (args_info->seed_given)
    write_into_file(outfile, "seed", args_info->seed_orig, 0);
  if (args_info->frac_error_bound_given)
    write_into_file(outfile, "frac-error-bound", args_info->frac_error_bound_orig, 0);
  if (args_info->surface_points_file_given)
    write_into_file(outfile, "surface-points-file", args_info->surface_points_file_orig, 0);
  if (args_info->interior_points_file_given)
    write_into_file(outfile, "interior-points-file", args_info->interior_points_file_orig, 0);
  if (args_info->print_counts_given)
    write_into_file(outfile, "print-counts", 0, 0 );
  if (args_info->print_benchmarks_given)
    write_into_file(outfile, "print-benchmarks", 0, 0 );
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  int result = EXIT_SUCCESS;

  if (cmdline_parser_required2(args_info, prog_name, 0) > 0)
    result = EXIT_FAILURE;

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required2 (struct gengetopt_args_info *args_info, const char *prog_name, const char *additional_error)
{
  int error_occurred = 0;
  FIX_UNUSED (additional_error);

  /* checks for required options */
  if (! args_info->input_file_given)
    {
      fprintf (stderr, "%s: '--input-file' ('-i') option required%s\n", prog_name, (additional_error ? additional_error : ""));
      error_occurred = 1;
    }
  
  
  /* checks for dependences among options */

  return error_occurred;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  char **string_field;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_DOUBLE:
    if (val) *((double *)field) = strtod (val, &stop_char);
    break;
  case ARG_LONGLONG:
#if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
    if (val) *((long long int*)field) = (long long int) strtoll (val, &stop_char, 0);
#else
    if (val) *((long *)field) = (long)strtol (val, &stop_char, 0);
#endif
    break;
  case ARG_STRING:
    if (val) {
      string_field = (char **)field;
      if (!no_free && *string_field)
        free (*string_field); /* free previous string */
      *string_field = gengetopt_strdup (val);
    }
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_DOUBLE:
  case ARG_LONGLONG:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error_occurred = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "full-help",	0, NULL, 0 },
        { "version",	0, NULL, 'V' },
        { "input-file",	1, NULL, 'i' },
        { "num-walks",	1, NULL, 0 },
        { "num-interior-samples",	1, NULL, 0 },
        { "max-rsd-capacitance",	1, NULL, 0 },
        { "max-rsd-polarizability",	1, NULL, 0 },
        { "max-rsd-volume",	1, NULL, 0 },
        { "min-num-walks",	1, NULL, 0 },
        { "min-num-interior-samples",	1, NULL, 0 },
        { "compute-form",	0, NULL, 0 },
        { "num-threads",	1, NULL, 0 },
        { "seed",	1, NULL, 0 },
        { "frac-error-bound",	1, NULL, 0 },
        { "surface-points-file",	1, NULL, 0 },
        { "interior-points-file",	1, NULL, 0 },
        { "print-counts",	0, NULL, 0 },
        { "print-benchmarks",	0, NULL, 0 },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVi:", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'i':	/* Input file name.  */
        
        
          if (update_arg( (void *)&(args_info->input_file_arg), 
               &(args_info->input_file_orig), &(args_info->input_file_given),
              &(local_args_info.input_file_given), optarg, 0, 0, ARG_STRING,
              check_ambiguity, override, 0, 0,
              "input-file", 'i',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          if (strcmp (long_options[option_index].name, "full-help") == 0) {
            cmdline_parser_print_full_help ();
            cmdline_parser_free (&local_args_info);
            exit (EXIT_SUCCESS);
          }

          /* Number of walk-on-spheres walks to perform.  */
          if (strcmp (long_options[option_index].name, "num-walks") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->num_walks_arg), 
                 &(args_info->num_walks_orig), &(args_info->num_walks_given),
                &(local_args_info.num_walks_given), optarg, 0, 0, ARG_LONGLONG,
                check_ambiguity, override, 0, 0,
                "num-walks", '-',
                additional_error))
              goto failure;
          
          }
          /* Number of interior samples to take.  */
          else if (strcmp (long_options[option_index].name, "num-interior-samples") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->num_interior_samples_arg), 
                 &(args_info->num_interior_samples_orig), &(args_info->num_interior_samples_given),
                &(local_args_info.num_interior_samples_given), optarg, 0, 0, ARG_LONGLONG,
                check_ambiguity, override, 0, 0,
                "num-interior-samples", '-',
                additional_error))
              goto failure;
          
          }
          /* Perform walk-on-spheres walks until the relative standard deviation of the capacitance drops below this value.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
          else if (strcmp (long_options[option_index].name, "max-rsd-capacitance") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->max_rsd_capacitance_arg), 
                 &(args_info->max_rsd_capacitance_orig), &(args_info->max_rsd_capacitance_given),
                &(local_args_info.max_rsd_capacitance_given), optarg, 0, 0, ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "max-rsd-capacitance", '-',
                additional_error))
              goto failure;
          
          }
          /* Perform walk-on-spheres walks until the relative standard deviation of the mean electric polarizability drops below this value.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
          else if (strcmp (long_options[option_index].name, "max-rsd-polarizability") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->max_rsd_polarizability_arg), 
                 &(args_info->max_rsd_polarizability_orig), &(args_info->max_rsd_polarizability_given),
                &(local_args_info.max_rsd_polarizability_given), optarg, 0, 0, ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "max-rsd-polarizability", '-',
                additional_error))
              goto failure;
          
          }
          /* Take interior samples until the relative standard deviation of volume drops below this value.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
          else if (strcmp (long_options[option_index].name, "max-rsd-volume") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->max_rsd_volume_arg), 
                 &(args_info->max_rsd_volume_orig), &(args_info->max_rsd_volume_given),
                &(local_args_info.max_rsd_volume_given), optarg, 0, 0, ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "max-rsd-volume", '-',
                additional_error))
              goto failure;
          
          }
          /* Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions.  */
          else if (strcmp (long_options[option_index].name, "min-num-walks") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->min_num_walks_arg), 
                 &(args_info->min_num_walks_orig), &(args_info->min_num_walks_given),
                &(local_args_info.min_num_walks_given), optarg, 0, "1000", ARG_LONGLONG,
                check_ambiguity, override, 0, 0,
                "min-num-walks", '-',
                additional_error))
              goto failure;
          
          }
          /* Minimum number of interior samples to take when using max-rsd stopping conditions.  */
          else if (strcmp (long_options[option_index].name, "min-num-interior-samples") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->min_num_interior_samples_arg), 
                 &(args_info->min_num_interior_samples_orig), &(args_info->min_num_interior_samples_given),
                &(local_args_info.min_num_interior_samples_given), optarg, 0, "10000", ARG_LONGLONG,
                check_ambiguity, override, 0, 0,
                "min-num-interior-samples", '-',
                additional_error))
              goto failure;
          
          }
          /* Compute form factor.  */
          else if (strcmp (long_options[option_index].name, "compute-form") == 0)
          {
          
          
            if (update_arg( 0 , 
                 0 , &(args_info->compute_form_given),
                &(local_args_info.compute_form_given), optarg, 0, 0, ARG_NO,
                check_ambiguity, override, 0, 0,
                "compute-form", '-',
                additional_error))
              goto failure;
          
          }
          /* Number of threads to use  (default=Number of logical cores).  */
          else if (strcmp (long_options[option_index].name, "num-threads") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->num_threads_arg), 
                 &(args_info->num_threads_orig), &(args_info->num_threads_given),
                &(local_args_info.num_threads_given), optarg, 0, 0, ARG_INT,
                check_ambiguity, override, 0, 0,
                "num-threads", '-',
                additional_error))
              goto failure;
          
          }
          /* Seed for the random number generator  (default=Randomly set).  */
          else if (strcmp (long_options[option_index].name, "seed") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->seed_arg), 
                 &(args_info->seed_orig), &(args_info->seed_given),
                &(local_args_info.seed_given), optarg, 0, 0, ARG_INT,
                check_ambiguity, override, 0, 0,
                "seed", '-',
                additional_error))
              goto failure;
          
          }
          /* Fractional error bound for nearest neighbor search.  */
          else if (strcmp (long_options[option_index].name, "frac-error-bound") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->frac_error_bound_arg), 
                 &(args_info->frac_error_bound_orig), &(args_info->frac_error_bound_given),
                &(local_args_info.frac_error_bound_given), optarg, 0, "0", ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "frac-error-bound", '-',
                additional_error))
              goto failure;
          
          }
          /* Name of file for writing the surface points from Walk-on-Spheres.  */
          else if (strcmp (long_options[option_index].name, "surface-points-file") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->surface_points_file_arg), 
                 &(args_info->surface_points_file_orig), &(args_info->surface_points_file_given),
                &(local_args_info.surface_points_file_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "surface-points-file", '-',
                additional_error))
              goto failure;
          
          }
          /* Name of file for writing the interior sample points.  */
          else if (strcmp (long_options[option_index].name, "interior-points-file") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->interior_points_file_arg), 
                 &(args_info->interior_points_file_orig), &(args_info->interior_points_file_given),
                &(local_args_info.interior_points_file_given), optarg, 0, 0, ARG_STRING,
                check_ambiguity, override, 0, 0,
                "interior-points-file", '-',
                additional_error))
              goto failure;
          
          }
          /* Print statistics related to counts of hit points.  */
          else if (strcmp (long_options[option_index].name, "print-counts") == 0)
          {
          
          
            if (update_arg( 0 , 
                 0 , &(args_info->print_counts_given),
                &(local_args_info.print_counts_given), optarg, 0, 0, ARG_NO,
                check_ambiguity, override, 0, 0,
                "print-counts", '-',
                additional_error))
              goto failure;
          
          }
          /* Print detailed RAM and timing information.  */
          else if (strcmp (long_options[option_index].name, "print-benchmarks") == 0)
          {
          
          
            if (update_arg( 0 , 
                 0 , &(args_info->print_benchmarks_given),
                &(local_args_info.print_benchmarks_given), optarg, 0, 0, ARG_NO,
                check_ambiguity, override, 0, 0,
                "print-benchmarks", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */



  if (check_required)
    {
      error_occurred += cmdline_parser_required2 (args_info, argv[0], additional_error);
    }

  cmdline_parser_release (&local_args_info);

  if ( error_occurred )
    return (EXIT_FAILURE);

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
