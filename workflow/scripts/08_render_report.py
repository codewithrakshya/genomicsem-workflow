#!/usr/bin/env python3
"""Build an offline HTML report from existing outputs (Python standard library only)."""
import argparse
import base64
import csv
import re
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from paths import results_dir


def read_rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def build_report(root, settings_path):
    settings = {r['key']: r['value'] for r in read_rows(root / settings_path)}
    traits_path = root / settings['traits_config']
    models_path = root / settings['models_config']
    traits, models = read_rows(traits_path), read_rows(models_path)
    analysis = settings['analysis_id']
    if Path(analysis).name != analysis or analysis in ('.', '..'):
        raise ValueError('analysis_id must be a filename component')
    result = results_dir(root)
    intermediate = root / '.work' / 'results' / result.name
    sections = []
    esc = lambda value: escape(str(value), quote=True)

    def table(path, limit=None):
        if not path.exists():
            return '<p class="missing">Not available: ' + esc(path.name) + '</p>'
        # Limit pilot previews without loading a potentially large result table.
        with path.open(newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader, [])
            rows = []
            total = 0
            for row in reader:
                total += 1
                if limit is None or len(rows) < limit:
                    rows.append(row)
        head = ''.join('<th scope="col">' + esc(c) + '</th>' for c in header)
        body = ''.join('<tr>' + ''.join('<td>' + esc(c) + '</td>' for c in row) + '</tr>' for row in rows)
        captions = {
            'genetic_covariances.tsv': 'Estimated genetic covariance matrix',
            'genetic_correlations.tsv': 'Estimated genetic correlations',
            'model_fit_comparison.tsv': 'Fit of the configured measurement models',
            'model_parameter_comparison.tsv': 'Factor loadings and residual variance estimates',
            'preparation_ldsc_qc.tsv': 'Variant preparation for LD score regression',
            'preparation_full_qc.tsv': 'Variant preparation for the full factor GWAS',
            'sampling_covariance_V.tsv': 'Sampling covariance of the genetic covariance estimates',
        }
        title = captions.get(path.name, path.stem.replace('_', ' '))
        if path.name.endswith('_factor_gwas_summary.tsv'):
            title = 'Factor-GWAS association and heterogeneity summary'
        elif path.name.endswith('_factor_gwas_pilot.tsv'):
            title = 'Pilot SNP-model results'
        caption = esc(title) + (f' — showing {len(rows)} of {total} rows' if limit else '')
        return f'<div class="table-scroll" tabindex="0"><table><caption>{caption}</caption><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'

    def text_file(path):
        if not path.exists():
            return '<p class="missing">Not available: ' + esc(path.name) + '</p>'
        text = path.read_text(errors='replace')
        return '<pre>' + esc(text or 'No warnings recorded in this file.') + '</pre>'

    def section(key, title, content):
        sections.append((key, title, f'<section id="{key}"><h2>{esc(title)}</h2>{content}</section>'))

    def fmt(value):
        try:
            return format(float(value), '.4g')
        except (ValueError, TypeError):
            return esc(value)

    primary = next(m for m in models if m['role'] == 'primary')
    section('overview', 'Introduction',
        '<p class="lead">Do the input phenotypes share a coherent genetic factor, and which common genetic variants are associated with that factor?</p>'
        '<p>This analysis combines GWAS summary statistics, not individual participant records. Genomic structural equation modeling (GenomicSEM) first models the genetic relationships among traits, then estimates SNP associations with the shared factor.</p>'
        f'<p>The analysis includes {len(traits)} traits and {len(models)} measurement models. The model labeled <strong>{esc(primary["label"])}</strong> is selected for the factor GWAS. “Primary” describes its pipeline role; it does not establish scientific validation.</p>'
        '<p><strong>Reading this report.</strong> Begin with the findings, then review the study definitions, genetic correlations, and model assumptions before interpreting SNP results. Expandable sections contain detailed estimates and model syntax. This report summarizes available outputs; it does not present independent replication.</p>')

    # Derive the findings from result tables so the report updates with each analysis.
    findings = []
    labels = {trait['trait']: trait.get('display_name', trait['trait']) for trait in traits}
    correlation_path = intermediate / 'genetic_correlations.tsv'
    if correlation_path.exists():
        correlations = read_rows(correlation_path)
        pairs = []
        for i, row in enumerate(correlations):
            for other in correlations[i + 1:]:
                value = row.get(other['trait'])
                if value and value != 'NA':
                    pairs.append(f"{esc(labels.get(row['trait'], row['trait']))}–{esc(labels.get(other['trait'], other['trait']))}: <strong>{float(value):.3f}</strong>")
        if pairs:
            findings.append('Estimated genetic correlations: ' + '; '.join(pairs) + '. These are point estimates of genetic sharing.')
    parameter_path = intermediate / 'model_parameter_comparison.tsv'
    if parameter_path.exists():
        parameters = read_rows(parameter_path)
        for model in models:
            column = 'Unstand_Est_' + model['label']
            negative = [row for row in parameters if row['op'] == '~~' and row['lhs'] == row['rhs']
                        and row.get(column) not in (None, '', 'NA') and float(row[column]) < 0]
            for row in negative:
                findings.append(f"The {esc(model['label'])} model estimates a negative variance for {esc(labels.get(row['lhs'], row['lhs']))} "
                                f"(<strong>{float(row[column]):.4f}</strong>), indicating an inadmissible solution. A boundary constraint requires scientific justification.")
    fit_path = intermediate / 'model_fit_comparison.tsv'
    if fit_path.exists():
        for fit in read_rows(fit_path):
            if fit['model'] == primary['label']:
                findings.append(f"Primary model ({esc(primary['label'])}): χ² = {fmt(fit.get('chisq', 'NA'))}, df = {fmt(fit.get('df', 'NA'))}, "
                                f"P = {fmt(fit.get('p_chisq', 'NA'))}; CFI = {fmt(fit.get('CFI', 'NA'))}; SRMR = {fmt(fit.get('SRMR', 'NA'))}. "
                                'Interpret these statistics together with the model constraints and residual estimates.')
    summary_path = result / f'{analysis}_factor_gwas_summary.tsv'
    if summary_path.exists():
        summary = {row['metric']: row['value'] for row in read_rows(summary_path)}
        def count(key):
            return format(int(summary[key]), ',') if key in summary else 'not available'
        findings.append(f"The factor GWAS analyzed <strong>{count('variants')}</strong> variants: "
                        f"<strong>{count('genome_wide_significant_rows')}</strong> association rows and "
                        f"<strong>{count('q_significant_rows')}</strong> heterogeneity rows met P &lt; 5 × 10⁻⁸. "
                        'Associated rows are not independent loci.')
    section('findings', 'Summary of findings', '<ul>' + ''.join('<li>' + item + '</li>' for item in findings) + '</ul>'
            if findings else '<p>Results are not yet available.</p>')


    stages = [
        ('00_prepare_sumstats.py', 'Python', 'Standardize columns, apply configured effect multipliers, remove invalid/duplicate variants, and select HapMap3 SNPs for LDSC.', 'preparation_ldsc_qc.tsv'),
        ('01_munge.R', 'GenomicSEM · munge()', 'Prepare and filter each trait’s summary statistics for multivariable LD score regression.', None),
        ('02_ldsc.R', 'GenomicSEM · ldsc()', 'Estimate genetic covariance S and sampling covariance V, accounting for correlated sampling error; derive genetic correlations.', 'genetic_covariances.tsv'),
        ('03_fit_model.R', 'GenomicSEM · usermodel()', 'Fit every configured measurement model using DWLS.', 'model_parameter_comparison.tsv'),
        ('04_compare_models.R', 'R', 'Collect fit statistics and parameters; flag negative residual variances.', 'model_fit_comparison.tsv'),
        ('05_visualize_results.R', 'R', 'Plot genetic correlations, factor loadings, and residual variances.', 'genomicsem_summary.png'),
        ('06_prepare_factor_gwas_pilot.py', 'Python', 'Select a deterministic pilot subset from the first prepared trait, then subset all traits and the reference.', None),
        ('07_factor_gwas.R · pilot', 'GenomicSEM · sumstats(), userGWAS()', 'Harmonize pilot SNPs; estimate factor associations and heterogeneity using iterative DWLS.', f'{analysis}_factor_gwas_pilot.tsv'),
        ('00_prepare_sumstats.py · full', 'Python', 'Prepare genome-wide inputs from raw files without the HapMap3 restriction.', 'preparation_full_qc.tsv'),
        ('07_factor_gwas.R · full', 'GenomicSEM · sumstats(), userGWAS()', 'Harmonize genome-wide SNP inputs; run analytic factor GWAS using the configured primary model.', f'{analysis}_factor_gwas_summary.tsv'),
    ]
    steps = '<ol class="steps">'
    for name, engine, description, artifact in stages:
        steps += f'<li><code>{esc(name)}</code><span class="engine">{esc(engine)}</span><p>{esc(description)}</p></li>'
    section('workflow', 'Statistical methods', steps + '</ol><p>The primary model is prespecified in the configuration; '
        'it is not automatically selected by fit statistics. The pilot checks SNP-model execution before the full factor GWAS.</p>')

    study_columns = [('display_name', 'Trait'), ('type', 'Phenotype type'), ('N', 'Configured sample size'),
                     ('ancestry', 'Ancestry'), ('build', 'Genome build'), ('effect_multiplier', 'Effect multiplier'),
                     ('sample_prev', 'Sample prevalence'), ('population_prev', 'Population prevalence')]
    study_table = '<div class="table-scroll"><table><caption>Study characteristics and modeling assumptions</caption><thead><tr>' + ''.join('<th scope="col">' + label + '</th>' for _, label in study_columns) + '</tr></thead><tbody>'
    study_table += ''.join('<tr>' + ''.join('<td>' + esc(trait.get(key, 'NA')) + '</td>' for key, _ in study_columns) + '</tr>' for trait in traits)
    study_table += '</tbody></table></div>'
    cohort_notes = ''
    cohort_path = root / settings['cohort_metadata']
    if cohort_path.exists():
        for row in read_rows(cohort_path):
            if row['trait'] in {t['trait'] for t in traits}:
                cohort_notes += '<p><strong>' + esc(row['study']) + '.</strong> ' + esc(row['phenotype']) + '. ' + esc(row['measurement']) + '; cohorts: ' + esc(row['cohorts']) + '. Reported sample: ' + esc(row['sample_size']) + '.</p>'
    cohort_notes += '<p><strong>Direction and scale.</strong> Reversing an effect changes its direction; it does not change the phenotype definition or convert binary effects into continuous measurement units. Configured effective sample sizes for case–control data can differ from participant counts. These remain distinct phenotypes after direction alignment.</p>'
    section('inputs', 'Study populations and phenotypes',
        cohort_notes +
        '<p>Effect direction is explicitly configured through <code>effect_multiplier</code>; it is not inferred from the GWAS filename. '
        'A multiplier of −1 reverses the effect sign without changing its SE or P value. Binary prevalence values are modeling assumptions.</p>'
        + study_table
        + '<p>Preparation does not carry an INFO column into its standardized output, so the configured INFO threshold alone '
        'does not establish filtering on original imputation quality. All traits must already use compatible genome builds and ancestry references.</p>')
    reference_rows = [
        ('HapMap3 SNP list', settings['hm3_reference'], 'Preparation and munge()',
         'Selects a common, well-characterized SNP set for LD score regression. It defines the LDSC SNP subset; the full factor GWAS is not restricted to this list.'),
        ('LD-score reference', settings['ld_reference'], 'ldsc(): ld',
         'Provides precomputed linkage-disequilibrium scores: how strongly each SNP tags surrounding variation. LDSC uses these scores to estimate genetic variances and cross-trait covariances from GWAS summary statistics.'),
        ('Regression-weight reference', settings['weight_reference'], 'ldsc(): wld',
         'Provides LD information used to weight the regression, accounting for correlated and differently variable association statistics. The same directory can supply both LD scores and regression weights.'),
        ('1000 Genomes SNP / allele / MAF reference', settings['gwas_reference'], 'sumstats() before userGWAS()',
         'Supplies reference alleles and allele frequencies for SNP harmonization and SNP variance calculations. Consistent allele coding ensures a SNP effect has the same direction across traits; reference frequencies support conversion between SNP effects and covariances.'),
    ]
    reference_html = '<p>GenomicSEM uses GWAS summary statistics rather than individual genotypes. External reference resources supply the LD and allele-frequency information needed to combine those statistics.</p>'
    reference_html += '<div class="table-scroll"><table><caption>Reference resources used in the analysis</caption><thead><tr><th>Resource</th><th>Configured file / directory</th><th>Used by</th><th>Purpose</th></tr></thead><tbody>'
    reference_html += ''.join('<tr>' + ''.join('<td style="white-space:normal;min-width:150px">' + esc(value) + '</td>' for value in row) + '</tr>' for row in reference_rows)
    reference_html += '</tbody></table></div>'
    reference_html += '<p>The example configuration points to European LD-score resources and the GenomicSEM 1000 Genomes Phase 3 reference; these resources must be downloaded separately. '
    reference_html += 'Configured trait ancestry: ' + esc(', '.join(sorted({t.get('ancestry', 'unspecified') for t in traits}))) + '; genome build: ' + esc(', '.join(sorted({t.get('build', 'unspecified') for t in traits}))) + '. '
    reference_html += 'Reference ancestry and variant coordinates must be compatible with the input GWAS. Changing the configured paths requires checking the replacement reference metadata.</p>'
    reference_html += '<p>HapMap3 is a SNP selection resource, LD scores describe correlations between variants, and the GWAS reference supplies alleles and frequencies; these resources serve different purposes. '
    reference_html += 'The reference filename’s MAF threshold is separate from the configured SNP-level MAF filter (' + esc(settings['maf_filter']) + ').</p>'
    reference_html += '<p>Methods and resource sources: <a href="https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models">GenomicSEM genome-wide models tutorial</a>; '
    reference_html += '<a href="https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS">GenomicSEM common-factor GWAS tutorial</a>; '
    reference_html += '<a href="https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v">reference downloads linked by GenomicSEM</a>.</p>'
    section('resources', 'Reference resources and why they are needed', reference_html)
    section('qc', 'Preparation quality control', table(intermediate / 'preparation_ldsc_qc.tsv')
        + '<p>“Kept” counts precede subsequent munging/harmonization. In LDSC mode, variants outside HapMap3 are skipped '
        'before duplicate and invalid-row counting; these columns are not an exhaustive breakdown of all exclusions.</p>'
        + table(intermediate / 'preparation_full_qc.tsv'))
    section('ldsc', 'Genetic covariance and correlation',
        '<p><strong>Question.</strong> How strongly are the genetic influences on the traits related?</p>'
        '<p><strong>Method.</strong> Multivariable LD score regression estimates genetic covariance S and sampling covariance V. V represents uncertainty and correlated sampling error, including that arising from overlapping participants.</p>'
        '<p>The diagonal of S contains SNP-heritability estimates on the analysis scales; off-diagonal entries are genetic covariances. '
        'Genetic correlations standardize these covariances: r<sub>g</sub>(i,j) = S<sub>ij</sub> / √(S<sub>ii</sub>S<sub>jj</sub>). '
        'The displayed correlation table contains point estimates, not confidence intervals.</p>'
        + table(intermediate / 'genetic_covariances.tsv') + table(intermediate / 'genetic_correlations.tsv')
        + '<p><strong>Interpretation.</strong> Positive genetic correlations indicate aligned genetic influences; they do not establish causation or make phenotypes interchangeable. Differences between pairs require assessment of uncertainty.</p>'
        + '<details><summary>Sampling covariance V</summary>' + table(root / '.work/ldsc/sampling_covariance_V.tsv') + '</details>')
    model_content = ('<p><strong>Question.</strong> Can a shared genetic factor describe the observed trait relationships?</p>'
        '<p><strong>Method.</strong> Diagonally weighted least squares (DWLS) fits the configured models, with factor variance fixed to one. Each trait’s genetic component is represented by a factor loading plus residual genetic variation.</p>'
        '<p>Loadings relate each trait to the shared genetic factor. Residual variances represent genetic variation '
        'not explained by that factor. Negative variances are inadmissible; fixing one to zero is a substantive model assumption. '
        'A model with zero degrees of freedom cannot provide an overall fit test. Consider chi-square alongside approximate fit indices.</p>')
    zero_residuals = []
    for model in models:
        syntax = (root / model['model_file']).read_text()
        for line in syntax.splitlines():
            match = re.fullmatch(r'\s*([A-Za-z_][A-Za-z0-9_.]*)\s*~~\s*0(?:\.0+)?\s*\*\s*\1\s*', line.split('#', 1)[0])
            if match and match.group(1) in labels:
                zero_residuals.append(f"{esc(labels[match.group(1)])} ({esc(model['label'])})")
    if zero_residuals:
        model_content += (
            '<div class="notice"><strong>Zero-residual constraint: rationale and references</strong>'
            '<p>The configured models fix residual genetic variance to zero for: ' + ', '.join(zero_residuals) + '. '
            'Within these models, the affected trait has no residual genetic variance outside the specified factor structure.</p>'
            '<p>Published precedent exists: Gard et al. (2021), Figure 4, report fixing the MDD residual variance to zero '
            'in both their one-factor and two-factor GenomicSEM models. '
            '<a href="#ref-gard">[1]</a> This is precedent for the modeling approach, not evidence that the residual '
            'variance of a trait in this analysis is truly zero.</p>'
            '<p>Variance constraints have consequences for estimation and model-fit inference '
            '<a href="#ref-savalei">[2]</a>; negative variance estimates should also be investigated as possible '
            'indicators of model misspecification <a href="#ref-kolenikov">[3]</a>. '
            'Removing a negative estimate alone does not validate the constrained model. '
            'The current primary model is selected in configuration, and its factor-GWAS findings are conditional '
            'on its constraints. Their use as primary assumptions requires scientific justification.</p></div>')
    for model in models:
        model_content += f'<details><summary>{esc(model["label"])} · {esc(model["role"])}</summary>'
        model_content += text_file(root / model['model_file'])
        model_content += '<p>SNP regression: <code>' + esc(model['snp_regression']) + '</code></p></details>'
        model_content += '<details><summary>Recorded warnings: ' + esc(model['label']) + '</summary>'
        model_content += text_file(root / '.work/models' / model['label'] / 'model_warnings.txt') + '</details>'
    section('models', 'Factor-model results', model_content + table(intermediate / 'model_fit_comparison.tsv')
        + text_file(intermediate / 'model_comparison_summary.txt')
        + '<details><summary>All parameter estimates</summary>' + table(intermediate / 'model_parameter_comparison.tsv') + '</details>')
    figure = result / 'genomicsem_summary.png'
    if figure.exists():
        encoded = base64.b64encode(figure.read_bytes()).decode('ascii')
        figure_html = f'<figure><img src="data:image/png;base64,{encoded}" alt="Genetic correlation heatmap, factor loadings with confidence intervals, and residual genetic variances">'
        figure_html += '<figcaption>Panels show genetic correlations, standardized loadings '
        'with approximate 95% intervals, and residual variances for the baseline and primary models.</figcaption></figure>'
    else:
        figure_html = '<p class="missing">Summary figure not available. Run the report stage after model fitting.</p>'
    section('figures', 'Visual summary of genetic sharing and model fit', figure_html)
    section('gwas', 'SNP associations with the common factor',
        '<p><strong>Method.</strong> GenomicSEM harmonizes SNP inputs with sumstats(), then fits the factor-on-SNP regression using userGWAS(). The pilot uses iterative DWLS; the full analysis uses analytic estimation with the primary measurement model.</p>'
        '<p>Factor associations test SNP effects on the latent factor scale, not directly in phenotype units. '
        'Q-SNP / omnibus heterogeneity tests departures from the effect pattern predicted by the common factor. '
        'Genome-wide significant rows are not independent loci; LD clumping or conditional analysis is needed.</p>'
        + table(result / f'{analysis}_factor_gwas_summary.tsv')
        + '<details><summary>Pilot preview (first 10 rows)</summary>' + table(result / f'{analysis}_factor_gwas_pilot.tsv', 10) + '</details>'
        + '<p>Optional outputs are shown if present; their presence is not a validation of the measurement model. '
        'The complete compressed GWAS remains a separate file to keep this report portable.</p>')
    section('overlap', 'Sample overlap and interpretation',
        '<p>Studies can share participants, so their estimates cannot be treated as independent. Multivariable LDSC estimates correlated sampling error, carried into model fitting through V. Accounting for dependence does not make overlapping studies independent replications.</p>'
        + '<details><summary>Study-specific overlap notes and sources</summary>' + text_file(root / settings['cohort_overlap']) + '</details>')
    section('interpretation', 'Discussion and limitations',
        '<p><strong>What is established here.</strong> Genetic correlations describe sharing among measured traits. Model comparisons show how the configured factor structures represent that covariance. Factor-GWAS associations are conditional on the selected measurement model.</p>'
        '<p><strong>What needs a scientific decision.</strong> A negative residual variance is inadmissible. Fixing it to zero changes the model assumption; it does not demonstrate absence of trait-specific genetic influences. Published precedent and favorable approximate fit indices alone do not validate that choice.</p>'
        '<p><strong>Before discovery claims.</strong> Review input coding, sample-size and prevalence assumptions, LDSC diagnostics, and the boundary constraint with the study team. Resolve associated SNP rows into independent loci and inspect heterogeneity. Biological annotation and independent replication are subsequent tasks, not findings supplied by this report.</p>')
    section('glossary', 'Terms used in this report',
        '<dl><dt>Genetic factor</dt><dd>A model-defined component of genetic variation shared across the traits.</dd>'
        '<dt>Loading</dt><dd>The relationship between a trait’s genetic component and the factor.</dd>'
        '<dt>Residual genetic variance</dt><dd>Genetic variation not explained by the specified factor structure; it is not simply environmental noise.</dd>'
        '<dt>CFI / SRMR</dt><dd>Comparative Fit Index / Standardized Root Mean Square Residual. Higher CFI and lower SRMR generally indicate closer fit; neither validates a biological interpretation.</dd>'
        '<dt>Q-SNP / omnibus heterogeneity</dt><dd>A test of departures from the cross-trait SNP-effect pattern predicted by the factor. A nonsignificant result does not prove an exact common-factor mechanism.</dd>'
        '<dt>SNP row versus independent locus</dt><dd>Nearby variants can be correlated through linkage disequilibrium (LD). Multiple associated rows may represent a single locus.</dd></dl>')
    if zero_residuals:
        section('references', 'References for variance constraints',
            '<ol>'
            '<li id="ref-gard">Gard AM, Ware EB, Hyde LW, Schmitz LL, Faul J, Mitchell C. (2021). '
            '<em>Phenotypic and genetic markers of psychopathology in a population-based sample of older adults.</em> '
            'Translational Psychiatry, 11, 239. Figure 4 reports the zero-residual constraint. '
            '<a href="https://doi.org/10.1038/s41398-021-01354-2">doi:10.1038/s41398-021-01354-2</a>.</li>'
            '<li id="ref-savalei">Savalei V, Kolenikov S. (2008). '
            '<em>Constrained versus unconstrained estimation in structural equation modeling.</em> '
            'Psychological Methods, 13(2), 150–170. '
            '<a href="https://doi.org/10.1037/1082-989X.13.2.150">doi:10.1037/1082-989X.13.2.150</a>.</li>'
            '<li id="ref-kolenikov">Kolenikov S, Bollen KA. (2012). '
            '<em>Testing Negative Error Variances: Is a Heywood Case a Symptom of Misspecification?</em> '
            'Sociological Methods &amp; Research, 41(1), 124–167. '
            '<a href="https://doi.org/10.1177/0049124112442138">doi:10.1177/0049124112442138</a>.</li>'
            '</ol>')
    nav = ''.join(f'<a href="#{key}">{esc(title)}</a>' for key, title, _ in sections)
    generated = datetime.now(timezone.utc).isoformat(timespec='seconds')
    css = """
:root{color-scheme:light;--link:#337ab7}*{box-sizing:border-box}body{margin:0;background:#fff;color:#292929;font:16px/1.75 Georgia,"Times New Roman",serif}header{max-width:1120px;margin:0 auto;padding:52px 30px 28px 250px;border-bottom:1px solid #ddd}h1{font:normal 36px/1.18 Georgia,serif;margin:0 0 18px;letter-spacing:-.5px}header p{font:13px/1.6 Arial,sans-serif;color:#666;margin:8px 0}.layout{display:grid;grid-template-columns:190px minmax(0,1fr);gap:30px;max-width:1120px;margin:auto;padding:30px}nav{position:sticky;top:24px;align-self:start;max-height:90vh;overflow:auto;font:12px/1.5 Arial,sans-serif;border:1px solid #ddd;border-radius:4px;padding:10px 0}nav:before{content:"Contents";display:block;padding:4px 14px 10px;font-weight:bold;color:#444}nav a{display:block;padding:6px 14px;color:#555;text-decoration:none}nav a:hover,nav a:focus{color:var(--link);background:#f5f5f5}nav a:first-of-type{color:var(--link)}main{min-width:0;counter-reset:section table figure}section{margin-bottom:38px;scroll-margin-top:24px}h2{font:normal 26px/1.3 Georgia,serif;border-bottom:1px solid #ddd;padding-bottom:8px;margin:0 0 18px;counter-increment:section}h2:before{content:counter(section)". ";color:#555}p{margin:12px 0}.lead{font-size:18px}.notice{border-left:3px solid #bbb;padding:2px 0 2px 16px;margin:22px 0}.missing{font-style:italic;color:#666}a{color:var(--link);text-decoration:none}a:hover{text-decoration:underline}li{margin:8px 0}.table-scroll{overflow:auto;margin:22px 0;max-height:600px}table{border-collapse:collapse;width:100%;font:12px/1.6 Arial,sans-serif;font-variant-numeric:tabular-nums;counter-increment:table}caption{text-align:left;font:13px/1.5 Georgia,serif;color:#444;padding-bottom:8px}caption:before{content:"Table " counter(table)". ";font-weight:bold}th,td{padding:8px 10px;text-align:left;vertical-align:top;white-space:nowrap;border-bottom:1px solid #ddd}th{border-top:2px solid #555;border-bottom:1px solid #555;background:white;position:sticky;top:0}tbody tr:nth-child(even){background:#fafafa}code,pre{font:12px/1.6 Menlo,Consolas,monospace}code{color:#444}pre{padding:12px;background:#f5f5f5;border:1px solid #e5e5e5;border-radius:3px;white-space:pre-wrap;overflow-wrap:anywhere}details{margin:12px 0;border:1px solid #ddd;border-radius:3px;padding:10px 12px}summary{cursor:pointer;font:13px/1.5 Arial,sans-serif;color:var(--link)}.steps{padding-left:22px}.steps p{margin:4px 0}.engine{font:11px/1.5 Arial,sans-serif;color:#777;display:block}figure{margin:24px 0;counter-increment:figure}img{width:100%;height:auto}figcaption{font:14px/1.6 Georgia,serif;margin-top:12px;color:#444}figcaption:before{content:"Figure " counter(figure)". ";font-weight:bold}dt{font-weight:bold;margin-top:14px}dd{margin:3px 0 14px}@media(max-width:800px){header{padding:28px 20px}.layout{display:block;padding:20px}nav{position:static;max-height:none;margin-bottom:28px;columns:2}h1{font-size:29px}h2{font-size:23px}}@media print{nav{display:none}header{padding:0 0 20px}.layout{display:block;padding:0}body{font-size:11pt}h2{break-after:avoid}.table-scroll{max-height:none;overflow:visible}th,td{white-space:normal;overflow-wrap:anywhere}figure{break-inside:avoid}a{color:black}details{border:0;padding:0}}
"""
    document = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
    document += f'<title>{esc(analysis)} · GenomicSEM report</title><style>{css}</style></head><body>'
    document += f'<header><h1>{esc(analysis.replace("_", " ").title())} — GenomicSEM analysis</h1><p>Genetic architecture, factor-model findings, and SNP associations · Updated {generated}</p></header>'
    document += f'<div class="layout"><nav aria-label="Table of contents">{nav}</nav><main>'
    document += ''.join(content for _, _, content in sections) + '</main></div></body></html>'
    result.mkdir(parents=True, exist_ok=True)
    destination = result / f'{analysis}_report.html'
    destination.write_text(document, encoding='utf-8')
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', default='.')
    parser.add_argument('--settings', default='.work/config/workflow.tsv')
    args = parser.parse_args()
    print('HTML report written to', build_report(Path(args.project_root).resolve(), Path(args.settings)))
