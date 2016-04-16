//
// Copyright (C) 2016 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

#include <catch.hpp>
#include <linalgwrap/Subscribable.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <memory>
#include <rapidcheck.h>
#include <rapidcheck/state.h>

// have an extra verbose output for rapidcheck function tests:
//#define HAVE_SUBSCRIPTION_RC_CLASSIFY

namespace linalgwrap {
namespace tests {
using namespace rc;

namespace subscription_tests {
/** An extremely simple class which is also
 *  subscribable */
struct SimpleSubscribable : public Subscribable {
    int data;
    SimpleSubscribable() : data(0){};
    explicit SimpleSubscribable(int d) : data{d} {};

    bool operator==(const SimpleSubscribable& other) const {
        return data == other.data;
    }
    bool operator!=(const SimpleSubscribable& other) const {
        return !operator==(other);
    }
};

//
// System under test (sut) and model
//

/** Note: The rc::state::check function makes vivid use of copying these
 *        classes around. So all state which is modified between runs
 *        should be living on the stack of this object */
template <typename SubscribableType>
struct SubscriptionSystemUnderTest {
    typedef SubscriptionPointer<SubscribableType> subscription_ptr_type;

    //
    // Members
    //
    /** The subscribable objects */
    std::vector<SubscribableType> objects;

    /** The subscription pointers */
    std::vector<subscription_ptr_type> pointers;

    //
    // Constructor and destructor
    //
    /** Generate a system under test with a few objects but no pointers. */
    SubscriptionSystemUnderTest(const std::vector<SubscribableType>& objects_)
          : objects{objects_}, pointers{} {}

    ~SubscriptionSystemUnderTest() {
        // Make sure the pointers get destroyed before the objects
        pointers.clear();

        // Now destroy the objects:
        objects.clear();
    }

    /** Make a full proper copy of the system including a proper rearrangement
     * of the pointers */
    SubscriptionSystemUnderTest(const SubscriptionSystemUnderTest& other)
          : objects(other.objects) {
        sensible_copy_pointers_from(other);
    }

    SubscriptionSystemUnderTest& operator=(
          const SubscriptionSystemUnderTest& other) {
        // Clear my content:
        pointers.clear();
        objects.clear();

        // Copy stuff:
        objects = other.objects;
        sensible_copy_pointers_from(other);
    }

    SubscriptionSystemUnderTest(SubscriptionSystemUnderTest&&) = delete;

    constexpr static size_t invalid_index = size_t(-1);
    size_t pointed_object_index_of(const subscription_ptr_type& ptr) const {
        for (size_t i = 0; i < objects.size(); ++i) {
            if (ptr.get() == &objects[i]) return i;
        }
        return invalid_index;
    }

  private:
    /** Copy the pointers from another object in a sensible way, i.e. take care
     * the actually point to the appropriate objects in this class.
     */
    void sensible_copy_pointers_from(const SubscriptionSystemUnderTest& other) {

        std::cout << "copy before" << other.pointers.size() << std::endl;

        pointers.reserve(other.pointers.size());
        for (auto itptr = std::begin(other.pointers);
             itptr != std::endl(other.pointers); ++itptr) {
            SubscribableType* raw = itptr->get();

            if (raw == nullptr) {
                pointers.emplace_back(itptr->subscriber_id());
                continue;
            }

            // find out where we point to:
            auto res = std::find_if(
                  std::begin(other.objects), std::end(other.objects),
                  [&](const SubscribableType& obj) { return &obj == raw; });

            if (res == std::end(other.objects)) {
                // something went wrong when finding the pointer.
                // Point to nullpointer to trigger an error.
                std::cerr << "Error finding an object referenced by a pointer."
                          << std::endl;
                pointers.emplace_back(itptr->subscriber_id());
                continue;
            }

            // find the index of the object:
            size_t index = res - std::begin(other.objects);

            // point to the correct one:
            pointers.emplace_back(itptr->subscriber_id(), objects[index]);
        }

        std::cout << "copy after" << pointers.size() << std::endl;
    }
};

/* Model for testing subscribable - subscription system
 *
 * The idea is to have a playground with a few objects
 * and subscription pointers and a function to check
 * consistency.
 *
 * Note: This model only works consistenly if each string
 *       used to identify a pointer is in fact unique!
 * **/
template <typename SubscribableType>
struct SubscriptionModel {
    //! struct mimicing a Subscribable
    struct subscribable_model {
        /** Collection of strings giving us the list
         * of subscribing objects */
        std::list<std::string> subscribers;

        /** The actual data object */
        SubscribableType data;
    };

    constexpr static size_t invalid_index = size_t(-1);

    //! struct mimicing a Subscription Pointer
    struct pointer_model {
        //! id of the pointer
        std::string id;

        //! Object we point to
        size_t object_index;

        pointer_model() = default;
        explicit pointer_model(std::string pid)
              : id(pid), object_index(invalid_index){};
    };

    //
    // Data members
    //
    //! The list of subscribers in the model
    std::vector<pointer_model> pointers;

    //! The list of objects in the test system:
    std::vector<subscribable_model> objects;

    //
    // Constructor and destructor
    //

    /** Construct a model situation from a set of objects
     *
     * In the resulting model no pointers are present.
     * */
    SubscriptionModel(const std::vector<SubscribableType>& objects_)
          : pointers{}, objects{objects_.size()} {

        auto construct_model = [](const SubscribableType& o) {
            subscribable_model m;
            m.data = o;
            return m;
        };

        std::transform(std::begin(objects_), std::end(objects_),
                       std::begin(objects), construct_model);
    }

    /** Destructor */
    ~SubscriptionModel() {
        // first release subscriptions
        pointers.clear();

        // then release objects
        objects.clear();
    }

    /** Proper copy constructor taking care of sanely copying the data in
     *  the pointer_model and not destroying the reference count of the
     *  objects.*/
    SubscriptionModel(const SubscriptionModel& other)
          : pointers{other.pointers}, objects{other.objects} {}

    //! Check if pointer id has been used by one of the pointers before
    bool is_pointer_id_used(const std::string& id) const {
        // predicate to check that the id of a pointer agrees
        auto check_id = [&](const pointer_model& m) { return m.id == id; };

        // See if any of the pointers has the same id.
        return std::any_of(std::begin(pointers), std::end(pointers), check_id);
    }

    /** Return a random string that is guaranteed to not be used by any of
     * the pointers in the model. */
    std::string unused_pointer_id() const {
        // static size_t i = 0;
        // return std::to_string(i++);
        return *gen::suchThat<std::string>([&](const std::string& s) {
            // Invert the result of the is_pointer_id_used function
            return !is_pointer_id_used(s);
        });
    }

    //
    // Assertive checks
    //

    /** Assert that an assertion is true for all objects in this model
     *  versus the system under test.
     *
     *  The ObjectAssertion is supposed to take an subscribable_model
     *  and a SubscribableType object and assert via RC_ASSERT that
     *  the condition is fulfilled for each such pair.
     *
     *  It is asserted that both sut as well as this model have the
     *  same number of objects.
     */
    template <typename ObjectAssertion>
    void assert_object_predicate_vs_sut(
          const SubscriptionSystemUnderTest<SubscribableType>& sut,
          ObjectAssertion oass) const {
        RC_ASSERT(objects.size() == sut.objects.size());

        auto itmodel = std::begin(objects);
        auto itsut = std::begin(sut.objects);
        for (; itmodel != std::end(objects); ++itmodel, ++itsut) {
            oass(*itmodel, *itsut);
        }
    }

    /** Assert that all model objects have all ids of the system under test
     *  No ordering is assumed
     */
    void assert_have_all_ids_of(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
#ifdef DEBUG
        auto assertion_on_pair = [](const subscribable_model& mod_obj,
                                    const SubscribableType& sut_obj) {
            // get the list of subscribers from the sut.
            auto sut_subscriptions = sut_obj.subscribers();

            // go through all ids:
            for (const auto& id : sut_subscriptions) {
                // is it contained in any of the subscribers of the
                // corresponding model object?
                auto res = std::any_of(
                      std::begin(mod_obj.subscribers),
                      std::end(mod_obj.subscribers),
                      [&](const std::string& mod_id) { return mod_id == id; });

                // assert that it is:
                RC_ASSERT(res);
            }

        };

        // Call the lambda defined above:
        assert_object_predicate_vs_sut(sut, assertion_on_pair);
#endif
    }

    /** Assert that all objects from the sut have all ids of the corresponding
     *  model object.
     *  No particular ordering is assumed
     */
    void assert_all_ids_in(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
#ifdef DEBUG
        auto assertion_on_pair = [](const subscribable_model& mod_obj,
                                    const SubscribableType& sut_obj) {
            // get the list of subscribers from the sut.
            auto sut_subscriptions = sut_obj.subscribers();

            // go through all ids:
            for (const auto& id : mod_obj.subscribers) {
                // is it contained in any of the subscribers of the
                // corresponding sut object?
                auto res = std::any_of(
                      std::begin(sut_subscriptions),
                      std::end(sut_subscriptions),
                      [&](const std::string& mod_id) { return mod_id == id; });

                // assert that it is:
                RC_ASSERT(res);
            }
        };

        // Call the lambda defined above:
        assert_object_predicate_vs_sut(sut, assertion_on_pair);
#endif
    }

    /**  Assert all objects from the sut have exactly the ids of the
     * corresponding
      *  model object.
      *  No particular ordering is assumed
      */
    void assert_ids_equal_to(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
#ifdef DEBUG
        // If all ids from this model are in the sut ...
        assert_all_ids_in(sut);

        // ... and the reverse is true, we are done.
        assert_have_all_ids_of(sut);
#endif
    }

    /** Assert that the pointer are equivalent. */
    static void assert_equal_pointer(
          const SubscriptionSystemUnderTest<SubscribableType>& sut,
          const SubscriptionPointer<SubscribableType>& sp,
          const pointer_model& mp) {

        // The index of the object we point to:
        size_t mto = mp.object_index;

        if (mto == invalid_index) {
            // Invalid index in mto, check that we point to nullptr:
            RC_ASSERT(sp.get() == nullptr);
        } else {
            size_t sutto = sut.pointed_object_index_of(sp);

            // Check that we point to the same thing:
            RC_ASSERT(sutto == mto);
        }
    }

    void assert_all_pointer_in(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
        auto& model = *this;

        // check for appropriate sizes of the containers:
        RC_ASSERT(model.pointers.size() <= sut.pointers.size());

        auto itsut = std::begin(sut.pointers);
        for (auto itmodel = std::begin(model.pointers);
             itmodel != std::end(model.pointers); ++itmodel, ++itsut) {
            assert_equal_pointer(sut, *itsut, *itmodel);
        }
    }

    void assert_have_all_pointer_of(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
        auto& model = *this;

        // check for appropriate sizes of the containers:
        RC_ASSERT(model.pointers.size() >= sut.pointers.size());

        auto itmodel = std::begin(model.pointers);
        for (auto itsut = std::begin(sut.pointers);
             itsut != std::end(sut.pointers); ++itsut, ++itmodel) {
            assert_equal_pointer(sut, *itsut, *itmodel);
        }
    }

    /** assert that all pointers are equal except a certain index
     *  (which is left unchecked).
     *  Asserts that the same *number* of pointers is present in both
     *  the model and the sut. */
    void assert_pointer_equal_to_except(
          const SubscriptionSystemUnderTest<SubscribableType>& sut,
          size_t idx) const {
        auto& model = *this;

        // check for equal sizes of the containers:
        RC_ASSERT(model.pointers.size() == sut.pointers.size());

        for (size_t i = 0; i < model.pointers.size(); ++i) {
            if (i == idx) continue;
            assert_equal_pointer(sut, sut.pointers[i], model.pointers[i]);
        }
    }

    void assert_pointer_equal_to(
          const SubscriptionSystemUnderTest<SubscribableType>& sut) const {
        auto& model = *this;

        // check for equal sizes of the containers:
        RC_ASSERT(model.pointers.size() == sut.pointers.size());

        // check that all our pointers are in the sut
        // together with the above check this assures, that
        // all pointers are equal
        assert_all_pointer_in(sut);
    }
};

//
// Operations and commands
//

/** Create a SubscriptionPointer with a random id pointing to a random
 * object */
template <typename SubscribableType>
struct CreateObjectPointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    std::string id;    //< id of the new pointer
    size_t obj_index;  //< object we point to
    CreateObjectPointer(const model_type& model)
          : id{model.unused_pointer_id()},
            // RHS of inRange is exclusive
            obj_index{*gen::inRange<size_t>(0, model.objects.size())} {};

    void apply(model_type& model) const override {
        // Make sure id has not been used yet
        RC_PRE(!model.is_pointer_id_used(id));

        // Create new model pointer
        typename model_type::pointer_model pm;
        pm.id = id;
        pm.object_index = obj_index;

        // Push id at the front:
        model.objects[obj_index].subscribers.push_front(id);

        // Push back a new model_pointer
        model.pointers.emplace_back(pm);
    }

    void run(const model_type& model, sut_type& sut) const override {
        // The model is the state before running the test
        // so this gives us the index, where the next pointer
        // will end up in both the model as well as the sut
        auto ptr_index = model.pointers.size();

#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "CreateObjectPointer");
#endif

        sut.pointers.push_back(
              make_subscription<SubscribableType>(sut.objects[obj_index], id));

        // check the id agrees.
        RC_ASSERT(sut.pointers[ptr_index].subscriber_id() == id);

        // Check that we point to the right thing
        RC_ASSERT(sut.pointers[ptr_index].get() == &sut.objects[obj_index]);

        // Check that the number of pointers has increased:
        RC_ASSERT(model.pointers.size() + 1 == sut.pointers.size());

        // Check that all the objects stiff have at least the references
        // we expect from the model.
        model.assert_all_ids_in(sut);

#ifdef DEBUG
        // Check that the object made note of us:
        auto obj_subscribers = sut.objects[obj_index].subscribers();
        auto res = std::any_of(
              std::begin(obj_subscribers), std::end(obj_subscribers),
              [&](const std::string& oid) { return oid == id; });
        RC_ASSERT(res);
#endif

        // check that the other pointers have not changed.
        model.assert_all_pointer_in(sut);
    }

    void show(std::ostream& os) const override {
        os << "Create SubscriptionPointer (id= " << id << ", to= " << obj_index
           << ")";
    }
};  // CreateObjectPointer

/** Create an empty SubscriptionPointer with a random id */
template <typename SubscribableType>
struct CreateEmptyPointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    std::string id;
    CreateEmptyPointer(const model_type& model)
          : id{model.unused_pointer_id()} {};

    void apply(model_type& model) const override {
        // Make sure id has not been used yet
        RC_PRE(!model.is_pointer_id_used(id));

        // Push back a new model_pointer
        model.pointers.emplace_back(id);
    }

    void run(const model_type& model, sut_type& sut) const override {
        // The model is the state before running the test
        // so this gives us the index, where the next pointer
        // will end up in both the model as well as the sut
        auto ptr_index = model.pointers.size();

#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "CreateEmptyPointer");
#endif

        sut.pointers.emplace_back(id);

        // check the id agrees.
        RC_ASSERT(sut.pointers[ptr_index].subscriber_id() == id);

        // Check that we point to null:
        RC_ASSERT(sut.pointers[ptr_index].get() == nullptr);

        // Check that the number of pointers has increased:
        RC_ASSERT(model.pointers.size() + 1 == sut.pointers.size());

        // Check that nothing has changed in the ids noted by the objects.
        model.assert_ids_equal_to(sut);

        // check that the other pointers have not changed.
        model.assert_all_pointer_in(sut);
    }

    void show(std::ostream& os) const override {
        os << "Create empty SubscriptionPointer (id=" << id << ", to=Null)";
    }
};  // CreateEmptyPointer

/** Take a random SubscriptionPointer and let it point somewhere else */
template <typename SubscribableType>
struct RedirectPointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    size_t ptr_index;  //< pointer we change
    size_t obj_index;  //< new object we point to
    RedirectPointer(const model_type& model)
          : ptr_index{*gen::inRange<size_t>(0, model.pointers.size())},
            obj_index{*gen::inRange<size_t>(0, model.objects.size())} {};

    void apply(model_type& model) const override {
        // check we have enough pointers to remove this one:
        RC_PRE(model.pointers.size() > ptr_index);
        RC_PRE(model.pointers[ptr_index].object_index != obj_index);

        typename model_type::pointer_model& pm = model.pointers[ptr_index];
        std::string id = pm.id;

        // delete the id from the currently pointed-to object subscribers list
        if (pm.object_index != model_type::invalid_index) {
            auto& subscr = model.objects[pm.object_index].subscribers;
            for (auto it = std::begin(subscr); it != std::end(subscr); ++it) {
                if (*it == id) {
                    subscr.erase(it);
                    break;
                }
            }
        }

        // Reset the pointer:
        pm.object_index = obj_index;

        // Push id at the front:
        model.objects[obj_index].subscribers.push_front(id);
    }

    void run(const model_type& model, sut_type& sut) const override {
        sut.pointers[ptr_index].reset(sut.objects[obj_index]);

#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "RedirectPointer");
#endif

        // Check that we point to the right thing
        RC_ASSERT(sut.pointers[ptr_index].get() == &sut.objects[obj_index]);

#ifdef DEBUG
        // Check that the object we pointed to does not
        // contain our id any more.
        // Only do this if we point to something
        auto& pm = model.pointers[ptr_index];
        if (pm.object_index != model_type::invalid_index) {
            auto subscribers = sut.objects[pm.object_index].subscribers();
            auto res = std::find(std::begin(subscribers), std::end(subscribers),
                                 pm.id);
            RC_ASSERT(res == std::end(subscribers));
        }

        // Check that the new object made note of us:
        auto obj_subscribers = sut.objects[obj_index].subscribers();
        auto res = std::any_of(
              std::begin(obj_subscribers), std::end(obj_subscribers),
              [&](const std::string& oid) { return oid == pm.id; });
        RC_ASSERT(res);

        // Check that the subscriptions agree
        auto assertion_on_pair = [&](
              const typename model_type::subscribable_model& mod_obj,
              const SubscribableType& sut_obj) {
            // get the list of subscribers from the sut.
            auto sut_subscriptions = sut_obj.subscribers();

            // go through all ids:
            for (const auto& id : mod_obj.subscribers) {
                // ignore the id we added fresh about:
                if (id == pm.id) continue;

                // is it contained in any of the subscribers of the
                // corresponding sut object?
                auto res = std::any_of(
                      std::begin(sut_subscriptions),
                      std::end(sut_subscriptions),
                      [&](const std::string& mod_id) { return mod_id == id; });

                // assert that it is:
                RC_ASSERT(res);
            }
        };

        // Call the lambda defined above:
        model.assert_object_predicate_vs_sut(sut, assertion_on_pair);
#endif

        // check that all pointers are equivalent except the one we
        // just changed
        model.assert_pointer_equal_to_except(sut, ptr_index);
    }

    void show(std::ostream& os) const override {
        os << "Reset SubscriptionPointer (ptr_index= " << ptr_index << ")";
    }
};  // RedirectPointer

/** Reset the pointing of a random SubscriptionPointer */
template <typename SubscribableType>
struct ResetPointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    size_t ptr_index;
    ResetPointer(const model_type& model)
          : ptr_index{*gen::inRange<size_t>(0, model.pointers.size())} {};

    void apply(model_type& model) const override {
        // check we have enough pointers to remove this one:
        RC_PRE(model.pointers.size() > ptr_index);

        typename model_type::pointer_model& pm = model.pointers[ptr_index];
        std::string id = pm.id;

        if (pm.object_index != model_type::invalid_index) {
            // delete the id from the object subscribers list
            auto& subscr = model.objects[pm.object_index].subscribers;
            for (auto it = std::begin(subscr); it != std::end(subscr); ++it) {
                if (*it == id) {
                    subscr.erase(it);
                    break;
                }
            }
        }

        // Reset the pointer:
        pm.object_index = model_type::invalid_index;
    }

    void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "ResetPointer");
#endif

        sut.pointers[ptr_index].reset();

        // Check that the object we pointed to does not
        // contain our id any more.
        // Only do this if we point to something
        auto& pm = model.pointers[ptr_index];
        if (pm.object_index != model_type::invalid_index) {
            auto subscribers = sut.objects[pm.object_index].subscribers();
            auto res = std::find(std::begin(subscribers), std::end(subscribers),
                                 pm.id);
            RC_ASSERT(res == std::end(subscribers));
        }

        // check that we point nowhere
        RC_ASSERT(sut.pointers[ptr_index].get() == nullptr);

        // Check that all ids of the objects in sut are also
        // referenced by the objects in the model.
        model.assert_have_all_ids_of(sut);

        // check that all pointers are equivalent except the one we
        // just changed
        model.assert_pointer_equal_to_except(sut, ptr_index);
    }

    void show(std::ostream& os) const override {
        os << "Reset SubscriptionPointer (ptr_index= " << ptr_index << ")";
    }
};  // ResetPointer

/** Remove a random SubscriptionPointer */
template <typename SubscribableType>
struct RemovePointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    size_t ptr_index;
    RemovePointer(const model_type& model)
          : ptr_index{*gen::inRange<size_t>(0, model.pointers.size())} {};

    void apply(model_type& model) const override {
        // check we have enough pointers to remove this one:
        RC_PRE(model.pointers.size() > ptr_index);

        typename model_type::pointer_model& pm = model.pointers[ptr_index];
        std::string id = pm.id;

        if (pm.object_index != model_type::invalid_index) {
            // delete the id from the object subscribers list
            auto& subscr = model.objects[pm.object_index].subscribers;
            for (auto it = std::begin(subscr); it != std::end(subscr); ++it) {
                if (*it == id) {
                    subscr.erase(it);
                    break;
                }
            }
        }

        // Delete the pointer:
        model.pointers.erase(ptr_index + std::begin(model.pointers));
    }

    void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "RemovePointer");
#endif

        sut.pointers.erase(ptr_index + std::begin(sut.pointers));

        // Check that the object we pointed to does not
        // contain our id any more.
        // Only do this if we point to something
        auto& pm = model.pointers[ptr_index];
        if (pm.object_index != model_type::invalid_index) {
            auto subscribers = sut.objects[pm.object_index].subscribers();
            auto res = std::find(std::begin(subscribers), std::end(subscribers),
                                 pm.id);
            RC_ASSERT(res == std::end(subscribers));
        }

        // Check that the number of pointers has decreased:
        RC_ASSERT(model.pointers.size() - 1 == sut.pointers.size());

        // Check that all ids of the objects in sut are also
        // referenced by the objects in the model.
        model.assert_have_all_ids_of(sut);

        // check that the other pointers have not changed.
        for (size_t i = 0; i < sut.pointers.size(); ++i) {
            auto i_mod = i;
            if (i >= ptr_index) i_mod = i + 1;

            model_type::assert_equal_pointer(sut, sut.pointers[i],
                                             model.pointers[i_mod]);
        }
    }

    void show(std::ostream& os) const override {
        os << "Delete SubscriptionPointer (ptr_index= " << ptr_index << ")";
    }
};  // RemovePointer

/** Copy a random SubscriptionPointer and delete the copy again.*/
template <typename SubscribableType>
struct CopyRemovePointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    size_t ptr_index;
    CopyRemovePointer(const model_type& model)
          : ptr_index{*gen::inRange<size_t>(0, model.pointers.size())} {};

    void apply(model_type& model) const override {
        // check we have enough pointers to copy this one:
        RC_PRE(model.pointers.size() > ptr_index);

        // Since our testing mechanism assumes, that now two pointers have the
        // same
        // id, we will copy the pointer, test the reference counting,
        // delete the copy and test the reference counting again.
        // So overall nothing should change.
    }

    void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "CopyRemovePointer");
#endif

        {
            // Make a copy:
            SubscriptionPointer<SubscribableType> copy(sut.pointers[ptr_index]);

            // check the id agrees.
            RC_ASSERT(sut.pointers[ptr_index].subscriber_id() ==
                      copy.subscriber_id());

            // Check that we point to the right thing
            RC_ASSERT(sut.pointers[ptr_index].get() == copy.get());

            // Check that all the objects still have at least the references
            // we expect from the model.
            model.assert_all_ids_in(sut);

#ifdef DEBUG
            if (model.pointers[ptr_index].object_index !=
                model_type::invalid_index) {
                // Extract object index:
                size_t obj_index =
                      sut.pointed_object_index_of(sut.pointers[ptr_index]);

                // Check that the object made note of us:
                auto obj_subscribers = sut.objects[obj_index].subscribers();
                auto res = std::count_if(std::begin(obj_subscribers),
                                         std::end(obj_subscribers),
                                         [&](const std::string& oid) {
                                             return oid == copy.subscriber_id();
                                         });
                RC_ASSERT(res == 2);
            }
#endif

            // check that the other pointers have not changed.
            model.assert_pointer_equal_to(sut);

        }  // Delete copy

        // check all objects now have exactly the references from
        // the model
        model.assert_ids_equal_to(sut);

        // check that the other pointers have not changed.
        model.assert_pointer_equal_to(sut);
    }

    void show(std::ostream& os) const override {
        os << "CopyRemove SubscriptionPointer (index= " << ptr_index << ")";
    }
};  // CopyRemovePointer

/** Assign a random SubscriptionPointer and delete the copy again */
template <typename SubscribableType>
struct AssignRemovePointer
      : rc::state::Command<SubscriptionModel<SubscribableType>,
                           SubscriptionSystemUnderTest<SubscribableType>> {
    typedef SubscriptionModel<SubscribableType> model_type;
    typedef SubscriptionSystemUnderTest<SubscribableType> sut_type;

    size_t ptr_index;
    AssignRemovePointer(const model_type& model)
          : ptr_index{*gen::inRange<size_t>(0, model.pointers.size())} {};

    void apply(model_type& model) const override {
        // check we have enough pointers to copy this one:
        RC_PRE(model.pointers.size() > ptr_index);

        // Since our testing mechanism assumes, that now two pointers have the
        // same id, we will copy-assign the pointer, test the reference
        // counting, delete the copy and test the reference counting again.
        // So overall nothing should change.
    }

    void run(const model_type& model, sut_type& sut) const override {
#ifdef HAVE_SUBSCRIPTION_RC_CLASSIFY
        RC_CLASSIFY(true, "AssignRemovePointer");
#endif

        {
            // Make an object and assign it:
            SubscriptionPointer<SubscribableType> copy("blubb");
            copy = sut.pointers[ptr_index];

            // check the id agrees.
            RC_ASSERT(sut.pointers[ptr_index].subscriber_id() ==
                      copy.subscriber_id());

            // Check that we point to the right thing
            RC_ASSERT(sut.pointers[ptr_index].get() == copy.get());

            // Check that all the objects still have at least the references
            // we expect from the model.
            model.assert_all_ids_in(sut);

#ifdef DEBUG
            if (model.pointers[ptr_index].object_index !=
                model_type::invalid_index) {
                // Extract object index:
                size_t obj_index =
                      sut.pointed_object_index_of(sut.pointers[ptr_index]);

                // Check that the object made note of us:
                auto obj_subscribers = sut.objects[obj_index].subscribers();
                auto res = std::count_if(std::begin(obj_subscribers),
                                         std::end(obj_subscribers),
                                         [&](const std::string& oid) {
                                             return oid == copy.subscriber_id();
                                         });
                RC_ASSERT(res == 2);
            }
#endif

            // check that the other pointers have not changed.
            model.assert_pointer_equal_to(sut);

        }  // Delete copy

        // check all objects now have exactly the references from
        // the model
        model.assert_ids_equal_to(sut);

        // check that the other pointers have not changed.
        model.assert_pointer_equal_to(sut);
    }

    void show(std::ostream& os) const override {
        os << "AssignRemove SubscriptionPointer (index= " << ptr_index << ")";
    }
};  // CopyRemovePointer

/** Move a random SubscriptionPointer to another, which is inserted into the
 * list of pointers and remove the old one */
// TODO implictly tested, but do explicit

/** Swap two random SubscriptionPointers and check that the reference counting
 * is the same, but their id and target have swapped */
// TODO implicitly checked via the assignment operator, but do explicit

}  // namespace subscription_tests

//
// ---------------------------------------------------------------
//

TEST_CASE("Subscription and SubscriptionPointer system", "[subscription]") {
    using namespace subscription_tests;

    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

    SECTION("SubscriptionPointer  gets the data from Subscribable") {
        SimpleSubscribable s{};
        s.data = 42;

        auto sptr = make_subscription(s, "Test");
        REQUIRE(sptr->data == s.data);
    }

    //
    // ---------------------------------------------------------------
    //

    SECTION("reset() and operator bool() of SubscriptionPointer.") {
        SimpleSubscribable s{};

        // create empty subscription
        SubscriptionPointer<SimpleSubscribable> sptr("Test");
        REQUIRE(!sptr);
        CHECK(s.n_subscriptions() == 0);

        // subscribe somewhere
        sptr.reset(s);
#ifdef DEBUG
        // n_subscriptions is always zero in RELEASE
        CHECK(s.n_subscriptions() == 1);
#endif

        REQUIRE(sptr);

        // reset subscription
        sptr.reset();
        REQUIRE(!sptr);
        CHECK(s.n_subscriptions() == 0);
    }

    //
    // ---------------------------------------------------------------
    //

    SECTION("Basic checks about Subscribables") {
        // create a Subscribable
        SimpleSubscribable s{};

        // change the data
        s.data = 42;

        // make two subscriptions:
        auto sub1 = make_subscription(s, "sub1");
        auto sub2 = make_subscription(s, "sub2");

#ifdef DEBUG
        SECTION("Test subscription was successful") {
            // Quick check for n_subscriptions
            REQUIRE(s.n_subscriptions() == 2);
            REQUIRE(s.subscribers()[0] == "sub2");
            REQUIRE(s.subscribers()[1] == "sub1");
        }
#endif

        SECTION("Test shared state of subscriptions") {
            // check if data is identical for starters.
            REQUIRE(sub2->data == sub1->data);

            // change data via subscription:
            sub1->data = 4;

            // check if the other got it as well
            REQUIRE(sub2->data == 4);
        }

        //
        // ---------------------------------------------------------------
        //

        SECTION("Copying subscribables") {
            // Copy s to s_copy
            SimpleSubscribable s_copy{s};

#ifdef DEBUG
            // Check number of subscriptions:
            CHECK(s.n_subscriptions() == 2);
            CHECK(s.subscribers()[0] == "sub2");
            CHECK(s.subscribers()[1] == "sub1");

            // Require the copy to be subscribed by
            // no one
            REQUIRE(s_copy.n_subscriptions() == 0);
#endif

            // Alter the data
            s.data = 43;

            // check consistency
            REQUIRE(sub1->data == 43);
            REQUIRE(s_copy.data == 42);
        }

        //
        // ---------------------------------------------------------------
        //

        SECTION("Assigning subscribables") {
            // Create an s_assigned
            SimpleSubscribable s_assigned;
            s_assigned.data = 1;

            // And subscribe to it once:
            auto sub3 = make_subscription(s_assigned, "sub3");

            // Assign to s_assigned:
            s_assigned = s;

            // Check assignment worked:
            REQUIRE(s_assigned.data == 42);

#ifdef DEBUG
            // Check number of subscriptions:
            CHECK(s.n_subscriptions() == 2);
            CHECK(s.subscribers()[0] == "sub2");
            CHECK(s.subscribers()[1] == "sub1");

            // Require the copy to be subscribed by
            // no one
            REQUIRE(s_assigned.n_subscriptions() == 1);
            REQUIRE(s_assigned.subscribers()[0] == "sub3");
#endif
            // Alter the data
            s_assigned.data = 5;

            // check consistency
            REQUIRE(sub1->data == 42);
            REQUIRE(sub3->data == 5);
            REQUIRE(s_assigned.data == 5);
        }
    }  // Basic checks about subscribables

    //
    // ---------------------------------------------------------------
    //

    SECTION("Random function test") {
        auto subscribing_to_different_objects_test =
              //[](std::vector<int> datas) {
              [] {
                  // Initial values:
                  std::vector<int> datas{1, 2, 3, 4, 5, 6, 7, 8, 9};

                  //
                  // The actual test
                  //
                  RC_PRE(datas.size() > size_t{0});

                  // The subscribable we will use for checking
                  typedef SimpleSubscribable subscribable_type;

                  // Allocate the test objects:
                  std::vector<subscribable_type> objects{datas.size()};
                  std::transform(std::begin(datas), std::end(datas),
                                 std::begin(objects),
                                 [](int d) { return SimpleSubscribable{d}; });

                  // Setup the model and initial system state
                  SubscriptionModel<subscribable_type> model{objects};
                  SubscriptionSystemUnderTest<subscribable_type> sut{objects};

                  // Typedef the operations
                  typedef CreateEmptyPointer<subscribable_type> op_CreateEmpty;
                  typedef CreateObjectPointer<subscribable_type> op_CreateObj;
                  typedef ResetPointer<subscribable_type> op_Reset;
                  typedef RedirectPointer<subscribable_type> op_Redirect;
                  typedef RemovePointer<subscribable_type> op_Remove;
                  typedef CopyRemovePointer<subscribable_type> op_CpRm;
                  typedef AssignRemovePointer<subscribable_type> op_AssignRm;

                  // Define generator for commands:
                  auto genCommands = state::gen::execOneOfWithArgs<
                        op_CreateObj, op_CreateEmpty, op_Reset, op_Remove,
                        op_Redirect, op_CpRm, op_AssignRm>;

                  // Run it through rapidcheck
                  state::check(model, sut, genCommands());
              };

        REQUIRE(rc::check(
              "Random function test of subscriber-subscription model.",
              subscribing_to_different_objects_test));
    }  // Random function test

    // TODO test comparison operators of SubscriptionPointer

}  // TEST_CASE
}  // namespace test
}  // namespace linalgwrap
