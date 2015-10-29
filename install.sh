# remember to run chmod 755 ./install.sh first time
# then run ./install.sh
echo ' ______________________________________________________ '
echo '|                                                      |'
echo '|                                                      |'
echo '|                                                      |'
echo '|                                                      |'
echo '|               installation begins  !                 |'
echo '|                                                      |'
echo '|                                                      |'
echo '|                                                      |'
echo '|______________________________________________________|'                                                
echo
echo 'Removing old buid'
rm -f ./transport
echo 'Starting build'
g++ -o transport src/main.cpp -std=c++11 || { 
echo "----- Build failed -----"
exit
}
echo 'Build complete'
echo 'Running transport code:'
echo '______________________________________________________'
echo 
./transport 