__author__ = 'Newton'
#
# Database Adapter for measurements
#
from contextlib import closing
import MySQLdb

host = "127.0.0.1"
user = "user"
password = "pass"
db = "alp"

test=True
def query_insert(query):
    if test:
        return
    try:
        con = MySQLdb.connect(host, user, password, db)
        with closing(con.cursor()) as cur:
            try:
                cur.execute(query)
                con.commit()
            except:
                con.rollback()
                import time
                time.sleep(5)
                query_insert(query)
                return
        con.close()
    except:
        print "Error with query:", query
        print "Could not read/write from database."

def query(sql):
    if test:
        return None
    try:
        con = MySQLdb.connect(host, user, password, db)
        with closing(con.cursor()) as cur:
            cur = con.cursor()
            cur.execute(sql)
            result = cur.fetchall()
        con.close()
    except:
        print "Error with query:", query
        print "Could not read/write from database."
        result = None

    return result
