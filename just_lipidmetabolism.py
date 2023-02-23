from oakvar import BasePostAggregator
from pathlib import Path
import sys
cur_path = str(Path(__file__).parent)
sys.path.append(cur_path)
import sqlite3
import lipidmetabolism_ref_homo


class CravatPostAggregator (BasePostAggregator):
    sql_insert = """ INSERT INTO lipid_metabolism (
                        rsid,
                        gene,
                        risk_allele,
                        genotype,
                        conclusion,
                        genotype_conclusion,
                        weight,
                        pmid,
                        population,
                        studydesign,
                        pvalue,
                        weightcolor
                    ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?) """
    ref_homo = lipidmetabolism_ref_homo.LipidRefHomo()

    def check(self):
        return True

    def setup (self):
        self.ref_homo.init(self, self.sql_insert)
        modules_path = str(Path(__file__).parent)
        sql_file = modules_path + "/data/lipid_metabolism.sqlite"
        if Path(sql_file).exists():
            self.lipid_conn = sqlite3.connect(sql_file)
            self.lipid_cursor = self.lipid_conn.cursor()

        self.result_path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.longevity_conn = sqlite3.connect(self.result_path)
        self.longevity_cursor = self.longevity_conn.cursor()
        sql_create = """ CREATE TABLE IF NOT EXISTS lipid_metabolism (
            id integer NOT NULL PRIMARY KEY,
            rsid text,
            gene text,
            risk_allele text,
            genotype text,
            conclusion text,
            genotype_conclusion text,
            weight float,
            pmid text,
            population text,
            studydesign text,
            pvalue text,
            weightcolor text
            )"""
        self.longevity_cursor.execute(sql_create)
        self.longevity_conn.commit()
        self.longevity_cursor.execute("DELETE FROM lipid_metabolism;")
        self.ref_homo.setup()


    def cleanup (self):
        if self.longevity_cursor is not None:
            self.longevity_cursor.close()
        if self.longevity_conn is not None:
            self.longevity_conn.commit()
            self.longevity_conn.close()
        if self.lipid_cursor is not None:
            self.lipid_cursor.close()
        if self.lipid_conn is not None:
            self.lipid_conn.close()
        return


    def get_color(self, w, scale = 1.5):
        w = float(w)
        if w < 0:
            w = w * -1
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = "ff" + color + color
        else:
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = color + "ff" + color

        return color


    def annotate (self, input_data):
        rsid = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        self.ref_homo.process_row(input_data)

        if not rsid.startswith('rs'):
            rsid = "rs" + rsid

        alt = input_data['base__alt_base']
        ref = input_data['base__ref_base']

        zygot = input_data['vcfinfo__zygosity']
        genome = alt + ref
        gen_set = {alt, ref}
        if zygot == 'hom':
            genome = alt + alt
            gen_set = {alt, alt}

        zygot:str = input_data['vcfinfo__zygosity']
        if zygot is None or zygot == "":
            zygot = "hom"

        query = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, rsid_conclusion, weight, " \
                " pmids, population, populations, p_value FROM rsids, studies, " \
                f" weight WHERE rsids.rsid = '{rsid}' AND weight.rsid = '{rsid}' AND studies.snp= '{rsid}' " \
                f" AND risk_allele='{alt}' AND zygot ='{zygot}' AND allele='{alt}'; "

        self.lipid_cursor.execute(query)
        rows = self.lipid_cursor.fetchall()

        if len(rows) == 0:
            return

        for row in rows:
            allele = row[0]
            row_gen = {row[2][0], row[2][1]}

            if gen_set == row_gen:
                task = (rsid, row[1], allele, genome, row[4], row[3], float(row[5]), row[6], row[7], row[8],
                        row[9], self.get_color(row[5], 0.6))
                self.longevity_cursor.execute(self.sql_insert, task)

        return {"col1":""}


    def postprocess(self):
        self.ref_homo.end()
        pass
