java -jar bin/ADB.jar \
-version_file version.xml \
-version_file /Users/jpilarski/intellij/feast/version.xml \
-version_file /Users/jpilarski/intellij/Beast2/version.xml \
-overwrite \
-seed 1 \
-loglevel debug \
-statefile example/inference.state \
example/inference.xml
