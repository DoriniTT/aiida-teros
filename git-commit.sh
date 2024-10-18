#git status
#echo 'What changes have you made?'
#read changes
#echo 'What files would you like to commit? (Type -A for all of them)'
#read files
#
#git add $files
#git commit -m "$changes"
#git pull 
#git push
#
git add -A
git commit -m "backup"
git pull
git push
