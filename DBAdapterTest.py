__author__ = 'Newton'
import DBAdapter as dba

result = dba.query("SELECT source FROM graphs WHERE source='aroadNetwork'")
print result
if result:
    print "we have a result!"
