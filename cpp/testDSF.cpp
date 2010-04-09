/*
 * testDSF.cpp
 *
 *   Created on: Mar 26, 2010
 *       Author: nikai
 *  Description: unit tests for DSF
 */

#include <iostream>
#include <boost/assign/std/list.hpp>
#include <boost/assign/std/set.hpp>
using namespace boost::assign;
#include <CppUnitLite/TestHarness.h>

#include "DSF.h"
#include "Key.h"

using namespace std;
using namespace gtsam;

/* ************************************************************************* */
TEST(DSF, makeSet) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	LONGS_EQUAL(1, dsf.size());
}

/* ************************************************************************* */
TEST(DSF, findSet) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	CHECK(dsf.findSet(5) != dsf.findSet(7));
}

/* ************************************************************************* */
TEST(DSF, makeUnion) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,7);
	CHECK(dsf.findSet(5) == dsf.findSet(7));
}

/* ************************************************************************* */
TEST(DSF, makeUnion2) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(7,5);
	CHECK(dsf.findSet(5) == dsf.findSet(7));
}

/* ************************************************************************* */
TEST(DSF, makeUnion3) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);
	dsf = dsf.makeUnion(6,7);
	CHECK(dsf.findSet(5) == dsf.findSet(7));
}

/* ************************************************************************* */
TEST(DSF, makePair) {
	DSFInt dsf;
	dsf = dsf.makePair(0, 1);
	dsf = dsf.makePair(1, 2);
	dsf = dsf.makePair(3, 2);
	CHECK(dsf.findSet(0) == dsf.findSet(3));
}

/* ************************************************************************* */
TEST(DSF, makeList) {
	DSFInt dsf;
	list<int> keys; keys += 5, 6, 7;
	dsf = dsf.makeList(keys);
	CHECK(dsf.findSet(5) == dsf.findSet(7));
}

/* ************************************************************************* */
TEST(DSF, numSets) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);
	LONGS_EQUAL(2, dsf.numSets());
}

/* ************************************************************************* */
TEST(DSF, sets) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeUnion(5,6);
	map<int, set<int> > sets = dsf.sets();
	LONGS_EQUAL(1, sets.size());

	set<int> expected; expected += 5, 6;
	CHECK(expected == sets[dsf.findSet(5)]);
}

/* ************************************************************************* */
TEST(DSF, sets2) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);
	dsf = dsf.makeUnion(6,7);
	map<int, set<int> > sets = dsf.sets();
	LONGS_EQUAL(1, sets.size());

	set<int> expected; expected += 5, 6, 7;
	CHECK(expected == sets[dsf.findSet(5)]);
}

/* ************************************************************************* */
TEST(DSF, sets3) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);
	map<int, set<int> > sets = dsf.sets();
	LONGS_EQUAL(2, sets.size());

	set<int> expected; expected += 5, 6;
	CHECK(expected == sets[dsf.findSet(5)]);
}

/* ************************************************************************* */
TEST(DSF, partition) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeUnion(5,6);

	list<int> keys; keys += 5;
	map<int, set<int> > partitions = dsf.partition(keys);
	LONGS_EQUAL(1, partitions.size());

	set<int> expected; expected += 5;
	CHECK(expected == partitions[dsf.findSet(5)]);
}

/* ************************************************************************* */
TEST(DSF, partition2) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);

	list<int> keys; keys += 7;
	map<int, set<int> > partitions = dsf.partition(keys);
	LONGS_EQUAL(1, partitions.size());

	set<int> expected; expected += 7;
	CHECK(expected == partitions[dsf.findSet(7)]);
}

/* ************************************************************************* */
TEST(DSF, partition3) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);

	list<int> keys; keys += 5, 7;
	map<int, set<int> > partitions = dsf.partition(keys);
	LONGS_EQUAL(2, partitions.size());

	set<int> expected; expected += 5;
	CHECK(expected == partitions[dsf.findSet(5)]);
}

/* ************************************************************************* */
TEST(DSF, set) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);
	set<int> set = dsf.set(5);
	LONGS_EQUAL(2, set.size());

	std::set<int> expected; expected += 5, 6;
	CHECK(expected == set);
}

/* ************************************************************************* */
int func(const int& a) { return a + 10; }
TEST(DSF, map) {
	DSFInt dsf;
	dsf = dsf.makeSet(5);
	dsf = dsf.makeSet(6);
	dsf = dsf.makeSet(7);
	dsf = dsf.makeUnion(5,6);

	DSFInt actual = dsf.map(&func);
	DSFInt expected;
	expected = expected.makeSet(15);
	expected = expected.makeSet(16);
	expected = expected.makeSet(17);
	expected = expected.makeUnion(15,16);
	CHECK(actual == expected);
}

/* ************************************************************************* */
TEST(DSF, flatten) {
	DSFInt dsf;
	dsf = dsf.makePair(1, 2);
	dsf = dsf.makePair(2, 3);
	dsf = dsf.makePair(5, 6);
	dsf = dsf.makePair(6, 7);
	dsf = dsf.makeUnion(2, 6);

	DSFInt actual = dsf.flatten();
	DSFInt expected;
	expected = expected.makePair(1, 2);
	expected = expected.makePair(1, 3);
	expected = expected.makePair(1, 5);
	expected = expected.makePair(1, 6);
	expected = expected.makePair(1, 7);
	CHECK(actual == expected);
}

/* ************************************************************************* */
TEST(DSF, flatten2) {
	Symbol x1('x',1);
	Symbol x2('x',2), x3('x',3), x4('x',4);
	list<Symbol> keys; keys += x1,x2,x3,x4;
	DSFSymbol dsf(keys);
	dsf = dsf.makeUnion(x1,x2);
	dsf = dsf.makeUnion(x3,x4);
	dsf = dsf.makeUnion(x1,x3);

	CHECK(dsf != dsf.flatten());

	DSFSymbol expected2;
	expected2 = expected2.makePair(x1, x2);
	expected2 = expected2.makePair(x1, x3);
	expected2 = expected2.makePair(x1, x4);
	CHECK(expected2 == dsf.flatten());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

