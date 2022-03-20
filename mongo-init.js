db.createUser(
    {
        user : "bot",
        pwd: "botpsswd",
        roles : [
            {
                "role" : "readWriteAnyDatabase",
                "db" : "admin"
            }
        ]
    }
);

 
print('===============JAVASCRIPT===============');
print('Count of rows in test collection: ' + db.test.count());
db.auth(
    {
        user : "bot",
        pwd: "botpsswd"
    }
);
print('==========bot authenticated successfully on admin db with readWriteAnyDatabase role==========');

db = db.getSiblingDB('test_db')

print('changed db to "test_db" (\'use "test_db"\' command)');

db.test_collection.insert({ myfield: 'test1', anotherfield: 'TEST1' });
db.test_collection.insert({ myfield: 'test2', anotherfield: 'TEST2' });

print('===============AFTER JS INSERT==========');
print('Count of rows in test collection: ' + db.test_collection.count());

alltest = db.test_collection.find();
while (alltest.hasNext()) {
  printjson(alltest.next());
}

if (db.test.count() == 2) print("==========all tests OK. proceeding further, mein fuhrer.==========");