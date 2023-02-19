import sqlite3
from sqlite3 import Error


class LipidRefHomo:

    rsid_map = {}

    def init(self, reporter, sql_insert):
        self.parent = reporter
        self.sql_insert = sql_insert


    def setup(self):
        sql = "SELECT rsid, risk_allele, weight FROM weight WHERE state = 'ref' AND zygosity = 'hom'"
        self.parent.lipid_cursor.execute(sql)
        rows = self.parent.lipid_cursor.fetchall()
        for row in rows:
            self.rsid_map[row[0]] = {'exist':True, 'allele':row[1], 'weight':row[2]}

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

    def process_row(self, row):
        rsid = str(row['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith('rs'):
            rsid = "rs"+rsid

        item = self.rsid_map.get(rsid)
        if item:
            self.rsid_map[rsid]['exist'] = False


    def end(self):
        for rsid in self.rsid_map:
            if self.rsid_map[rsid]['exist']:
                allele = self.rsid_map[rsid]['allele']
                genotype = allele+allele

                query = "SELECT rsids.risk_allele, gene, genotype, genotype_specific_conclusion, rsid_conclusion, weight, " \
                " pmids, population, populations, p_value FROM rsids, studies, " \
                f" weight WHERE rsids.rsid = '{rsid}' AND weight.rsid = '{rsid}' AND studies.snp= '{rsid}' AND genotype = '{genotype}';"

                self.parent.lipid_cursor.execute(query)
                row = self.parent.lipid_cursor.fetchone()
                if row:
                    task = (rsid, row[1], allele, genotype, row[4], row[3], float(row[5]), row[6], row[7], row[8],
                        row[9], self.get_color(row[5], 0.6))
                    self.parent.longevity_cursor.execute(self.sql_insert, task)