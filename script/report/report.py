import jinja2
import fire
from pathlib import Path, PurePath
import os
import glob
import sys

script_dir, _ = os.path.split(os.path.abspath(__file__))

env = jinja2.Environment(loader=jinja2.FileSystemLoader(
    searchpath='{}/template'.format(script_dir)
)
)
template = env.get_template('index.html')
# code that fills in display_dictionary with the values to send to the template


def table2dict(table_file, name, sep='\t'):
    table_dict = dict()
    with open(table_file) as table_inf:
        for n, eachline in enumerate(table_inf):
            eachline_inf = eachline.split(sep)
            if n == 0:
                label = '{}_header'.format(name)
                table_dict[label] = eachline_inf
            else:
                label = '{}_body'.format(name)
                table_dict.setdefault(label, []).append(eachline_inf)
    return table_dict


def plot2report(plot_path, outpath, plot_name=None):
    plots = glob.glob(str(plot_path))
    outpath = PurePath(outpath)
    if plots:
        plot = plots[0]
        plot_path = PurePath(plot)
        if plot_name is None:
            plot_name = plot_path.stem
        outfile_path = outpath / f'{plot_name}{plot_path.suffix}'
        os.system(f'cp {plot_path} {outfile_path}')
    else:
        sys.exit(f'Can not find file: {plot_path}')


def exom_report(result_dir, proj_name, report_dir=None):
    result_dir = Path(result_dir)
    if report_dir is None:
        report_dir = result_dir / 'report'
    else:
        report_dir = Path(report_dir)
    display_dictionary = {}
    display_dictionary['project_name'] = proj_name

    # add fastqc table
    qc_table = result_dir / 'stats/fastqc.summary.xls'
    display_dictionary.update(
        table2dict(qc_table, 'seq'))

    # add aligment table
    align_table = result_dir / 'stats/all_sample.mapping.xls'
    display_dictionary.update(
        table2dict(align_table, 'align'))

    # snp stats
    # summary
    snp_summary_table = result_dir / 'stats/overall.varSummary.txt'
    display_dictionary.update(
        table2dict(snp_summary_table, 'snp_summary'))

    snp_number_table = result_dir / 'stats/overall.varNum.txt'
    display_dictionary.update(
        table2dict(snp_number_table, 'snp_number'))

    snp_impact_table = result_dir / 'stats/overall.varImpact.txt'
    display_dictionary.update(
        table2dict(snp_impact_table, 'snp_impact'))

    snp_effect_table = result_dir / 'stats/overall.varEffects.txt'
    display_dictionary.update(
        table2dict(snp_effect_table, 'snp_effect'))

    snp_region_table = result_dir / 'stats/overall.varRegion.txt'
    display_dictionary.update(
        table2dict(snp_region_table, 'snp_region'))

    report_dir.mkdir(parents=True, exist_ok=True)
    os.system('cp -r {script_dir}/template/* {report_dir}'.format(
        script_dir=script_dir,
        report_dir=report_dir
    ))

    display_html = template.render(display_dictionary)
    report_html = report_dir / 'index.html'
    with open(report_html, 'w') as out_inf:
        out_inf.write(display_html)

    # plots
    report_plot_path = report_dir / 'imgs'
    mapping_plot = result_dir / 'stats/plot/mapping/Mapping_stats.png'
    plot2report(mapping_plot, report_plot_path)

    genome_cov_plot = result_dir / 'stats/plot/coverage/Reads_coverage_genome.png'
    plot2report(genome_cov_plot, report_plot_path)
    cds_cov_plot = result_dir / 'stats/plot/coverage/Reads_coverage_cds.png'
    plot2report(cds_cov_plot, report_plot_path)

    variant_summary_plot = result_dir / \
        'stats/plot/variants/Variant_stats_summary.png'
    plot2report(variant_summary_plot, report_plot_path)

    varType_plot = result_dir / 'stats/plot/variants/*varType.png'
    plot2report(varType_plot, report_plot_path, 'varType')
    varRegion_plot = result_dir / 'stats/plot/variants/*varRegion.png'
    plot2report(varRegion_plot, report_plot_path, 'varRegion')
    varEffects_high_plot = result_dir / 'stats/plot/variants/*varEffects-HIGH.png'
    plot2report(varEffects_high_plot, report_plot_path, 'varEffects-HIGH')
    varEffects_moderate_plot = result_dir / \
        'stats/plot/variants/*varEffects-MODERATE.png'
    plot2report(varEffects_moderate_plot,
                report_plot_path, 'varEffects-MODERATE')
    varEffects_low_plot = result_dir / 'stats/plot/variants/*varEffects-LOW.png'
    plot2report(varEffects_low_plot, report_plot_path, 'varEffects-LOW')
    varEffects_modifier_plot = result_dir / \
        'stats/plot/variants/*varEffects-MODIFIER.png'
    plot2report(varEffects_modifier_plot,
                report_plot_path, 'varEffects-MODIFIER')
    varImpact_plot = result_dir / \
        'stats/plot/variants/*varImpact.png'
    plot2report(varImpact_plot,
                report_plot_path, 'varImpact')

    deltaSNP_plot = result_dir / 'mapping/*deltaSNP.png'
    Gprime_plot = result_dir / 'mapping/*Gprime.png'
    negLog10Pval_plot = result_dir / 'mapping/*negLog10Pval.png'
    plot2report(deltaSNP_plot, report_plot_path, 'deltaSNP')
    plot2report(Gprime_plot, report_plot_path, 'Gprime')
    plot2report(negLog10Pval_plot, report_plot_path, 'negLog10Pval')


if __name__ == '__main__':
    fire.Fire(exom_report)
